#include <Eigen/Eigen>

#include "use-ikfom.hpp"
#include <typeinfo>

/*============================================定义一些地理坐标常数==============================================*/
#define OMEGA 7.2921151467e-5  // rotational angular velocity of the earth
#define C_SPEED 2.99792458e8

/* Physical parameters of the Earth, Sun and Moon  */
#define R_WGS84 6378137.0           /* Radius Earth [m]; WGS-84  */
#define B_WGS84 6356752.3141        /* Radius(B) Earth [m]; WGS-84  */
#define F_WGS84 1.0 / 298.257223563 /* Flattening; WGS-84   */
#define Omega_WGS 7.2921151467e-5   /*[rad/s], the earth rotation rate */
#define GM_Earth 398600.5e+9        /* [m^3/s^2]; WGS-84 */
#define GM_JGM3 398600.4418e+9      /* [m^3/s^2]; JGM3  */
#define E2_WGS84 0.0066943799901413
#define GammaA 9.7803253359
#define GammaB 9.8321849379

/* Physical parameters of the Earth, Sun and Moon  */
#define R_CGS2K 6378137.0           /* Radius Earth [m]; CGCS2000  */
#define F_CGS2K 1.0 / 298.257222101 /* Flattening; CGCS2000   */
#define E2_CGS2K 0.0066943800229008
#define GM_BDS 398600.4418e+9 /* [m^3/s^2]; CGCS2000  */
/*===========================================================================================================*/

/*******************************************************
 * @brief       计算卯酉圈曲率半径
 * @author      Sihan Wang
 * @paras[in]   fai（纬度，rad）
 * @paras[out]  none
 * @return      RN
 * @Date        20/5/2022
 *******************************************************/
double cal_RN(double fai) {
  double W = sqrt(1.0 - E2_WGS84 * sin(fai) * sin(fai));
  double RN = R_WGS84 / W;
  return RN;
}

/*******************************************************
 * @brief       计算子午圈曲率半径
 * @author      Sihan Wang
 * @paras[in]   fai（纬度，rad）
 * @paras[out]  none
 * @return      RM
 * @Date        20/5/2022
 *******************************************************/
double cal_RM(double fai) {
  double W = sqrt(1.0 - E2_WGS84 * sin(fai) * sin(fai));
  double RM = R_WGS84 * (1.0 - E2_WGS84) / (W * W * W);
  return RM;
}

/*******************************************************
 * @brief       反对称矩阵
 * @author      Sihan Wang
 * @paras[in]   三维向量Vec
 * @paras[out]  none
 * @return      反对称矩阵Skew
 * @Date        20/5/2022
 *******************************************************/
Matrix3d SkewMat(V3D Vec) {
  Matrix3d Skew = Matrix3d::Zero(3, 3);
  Skew(0, 1) = -Vec(2);
  Skew(0, 2) = Vec(1);
  Skew(1, 0) = Vec(2);
  Skew(1, 2) = -Vec(0);
  Skew(2, 0) = -Vec(1);
  Skew(2, 1) = Vec(0);
  return Skew;
}

class GNSSProcess {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    VD(6) Z;  // Z: 观测误差向量
    MD(6,6) R; // R: 观测值的方差矩阵
    MD(6,23) H;    // H: 观测值的系数矩阵

    GNSSProcess();
    ~GNSSProcess();

    void Process(const MeasureGroup &meas, esekfom::esekf<state_ikfom, 12, input_ikfom> &kf_state, shared_ptr<V3D> leverarm, const V3D &gyro);



  private:
    void GNSS_Init(const V3D &lla_cord);
    V3D ConvertToNED(const V3D &lla_cord);
    void CorrectLeverArm(const V3D &leverarm, const V3D &w_bi_b);

    /* 姿态变量 */
    bool gnss_need_init_ = true;
    V3D* init_pos;
    V3D imu_pos;
    V3D imu_vel;
    /* 定义地理坐标变量 */
    double lat, lon, alt;
    double vN, vE, vD;
    double RM, RN;
    V3D w_ei_n;  // 地球自转角速度在n系下的投影
    V3D w_ne_n;  // n系相对于e系的角速度在n系下的投影（牵连角速度）
    V3D w_ni_n;  // n系相对于i系的角速度在n系下的投影（牵连角速度）

    /* 定义姿态变量 */
    M3D C_b2n;
    M3D I33 = Matrix3d::Identity();
    MD(23,23) IFF = MatrixXd::Identity(23,23);
    /* 定义GNSS观测状态变量 */
    V3D gnss_lla_cord;
    V3D gnss_lla_var;
    V3D gnss_vel;
    V3D gnss_vel_var;
    V3D gnss_pos;
};

GNSSProcess::GNSSProcess()
    : gnss_need_init_(true), init_pos(new V3D)
{}

GNSSProcess::~GNSSProcess() {}

V3D GNSSProcess::ConvertToNED(const V3D &lla_cord) {
  V3D ned;
  ned.setZero();
  if (gnss_need_init_) {
    ROS_WARN("GNSS Convert To NED Cord Function Uninited!");
    return ned;
  }
  else if((lla_cord(0) < -90.0) || (lla_cord(0) > +90.0) || (lla_cord(1) < -180.0) || (lla_cord(1) > +180.0)) // 经纬度有效性检查
  {
    ROS_WARN("WGS lat or WGS lon out of range");
    return ned;
  }

  /* LLA 转 ECEF 坐标系 */
  double sin_lat = sin(deg2rad(lla_cord(0)));
  double cos_lat = cos(deg2rad(lla_cord(0)));
  double sin_lon = sin(deg2rad(lla_cord(1)));
  double cos_lon = cos(deg2rad(lla_cord(1)));
  double r_n = R_WGS84/sqrt(1.0 - E2_WGS84 * sin_lat * sin_lat);
  V3D xyz = {
    (r_n + lla_cord(2)) * cos_lat * cos_lon,
    (r_n + lla_cord(2)) * cos_lat * sin_lon,
    (r_n * (1.0 - E2_WGS84) + lla_cord(2)) * sin_lat
  };

  /* ECEF 转 NED 坐标系 */
  double cos_lat_ref = cos(deg2rad((*init_pos)(0)));
  double sin_lat_ref = sin(deg2rad((*init_pos)(0)));
  double cos_lon_ref = cos(deg2rad((*init_pos)(1)));
  double sin_lon_ref = sin(deg2rad((*init_pos)(1)));
  double r_n_ref = R_WGS84 / sqrt(1.0 - E2_WGS84 * sin_lat_ref * sin_lat_ref);
  // 先把初始参考坐标转换到ECEF坐标系下
  V3D ref_xyz = {
    (r_n_ref + (*init_pos)(2)) * cos_lat_ref * cos_lon_ref,
    (r_n_ref + (*init_pos)(2)) * cos_lat_ref * sin_lon_ref,
    (r_n_ref * (1.0 - E2_WGS84) + (*init_pos)(2)) * sin_lat_ref
  };

  V3D diff_xyz = xyz - ref_xyz;
  ned = {
    -sin_lat_ref * cos_lon_ref * diff_xyz(0) - sin_lat_ref * sin_lon_ref * diff_xyz(1) + cos_lat_ref * diff_xyz(2),
    cos_lat_ref * sin_lon_ref * diff_xyz(0) - cos_lat_ref * cos_lon_ref * diff_xyz(1),
    -sin_lon_ref * diff_xyz(0) + cos_lon_ref * diff_xyz(1)
  };
  return ned;
}

void GNSSProcess::GNSS_Init(const V3D &lla_cord) {
  if(lla_cord.sum() <= 1e-3)
  {
    ROS_WARN("GPS NavSat Unfixed!");
  }
  else
  {
    *init_pos = lla_cord;
    gnss_need_init_ = false;
    ROS_INFO("GNSS Initial Done");
  }
}
void GNSSProcess::CorrectLeverArm(const V3D &leverarm, const V3D &w_bi_b)
{
  M3D mat_n2l = Matrix3d::Zero(3, 3); //将北东地转化为LLH系
	mat_n2l(0, 0) = 1.0 / (RM + alt);
	mat_n2l(1, 1) = 1.0 / ((RN + alt) * cos(lat));
	mat_n2l(2, 2) = -1.0;

  M3D OMIGA_ni_n = SkewMat(w_ni_n);
	M3D OMIGA_bi_b = SkewMat(w_bi_b);
	/*位置改正*/
	gnss_pos = gnss_pos - mat_n2l * C_b2n * leverarm;
	/*速度改正*/
	gnss_vel = gnss_vel + C_b2n * OMIGA_bi_b * leverarm + OMIGA_ni_n * C_b2n * leverarm;
}

void GNSSProcess::Process(const MeasureGroup &meas, esekfom::esekf<state_ikfom, 12, input_ikfom> &kf_state, shared_ptr<V3D> leverarm, const V3D &gyro) {
  /* 获取当前系统状态变量 */
  state_ikfom init_state = kf_state.get_x(); //rot pos vel bias_gyro bias_acc gravity R_L p_L
  /* 获取当前系统协方差矩阵 */
  esekfom::esekf<state_ikfom, 12, input_ikfom>::cov init_P = kf_state.get_P(); // 24x24维

  imu_pos = init_state.pos;
  imu_vel = init_state.vel;
  // imu_eul = SO3ToEuler(init_state.rot);
  /* 获取当前系统状态变量 */
  uint8_t N = 1;
  for (const auto &gnss : meas.gnss)
  {
    const auto &lla_cord = gnss->segment<3>(1);
    const auto &lla_var = gnss->segment<3>(4);
    const auto &vel = gnss->segment<3>(7);
    const auto &vel_var = gnss->segment<3>(10);
    if(lla_cord.sum() <= 1e-3)
    {
      break;
    }
    gnss_lla_cord += (lla_cord - gnss_lla_cord) / N;
    gnss_lla_var  += (lla_var - gnss_lla_var) / N;
    gnss_vel      += (vel - gnss_vel) / N;
    gnss_vel_var  += (vel_var - gnss_vel_var) / N;
    // cout<<"acc norm: "<<cur_acc.norm()<<" "<<mean_acc.norm()<<endl;
    N++;
  }
  /* 若GNSS需要初始化 */
  if(gnss_need_init_) GNSS_Init(gnss_lla_cord);
  else
  {
    lat = imu_pos(0);
    lon = imu_pos(1);
    alt = imu_pos(2);

    vN = imu_vel(0);
    vE = imu_vel(1);
    vD = imu_vel(2);

    RM = cal_RM(lat);
    RN = cal_RN(lat);

    w_ei_n << Omega_WGS * cos(lat), 0, -1.0 * Omega_WGS * sin(lat);
    w_ne_n << vE / (RN + alt), -1.0 * vN / (RM + alt), -1.0 * vE * tan(lat) / (RN + alt);
    w_ni_n = w_ei_n + w_ne_n;

    C_b2n = init_state.rot.toRotationMatrix();
    gnss_pos = ConvertToNED(gnss_lla_cord);

    ROS_INFO("Lat:%0.6f Lon:%0.6f Alt:%0.6f", gnss_lla_cord(0), gnss_lla_cord(1), gnss_lla_cord(2));
    ROS_INFO("N:%0.6f E:%0.6f D:%0.6f", gnss_pos(0), gnss_pos(1), gnss_pos(2));

    CorrectLeverArm(*leverarm, gyro); //杆臂效应改正

    VD(23) x;
    x.setZero();
    x.block<3, 1>(0, 3) = init_state.pos;
    x.block<3, 1>(0, 6) = init_state.vel;

    // DR_inv 将 n 系下的北向、东向和垂向位置差异（单位 m）转化为纬度、经度和高程分量的差异
    M3D DR;
    DR.setZero();
    DR(0, 0) = RM + alt;
    DR(1, 1) = (RN + alt) * cos(lat);
    DR(2, 2) = -1;
    M3D DR_inv = DR.inverse();

    /***********************观测方程***************************/
    Z.setZero();
    Z.block<3, 1>(0, 0) = DR * (imu_pos - gnss_pos) + C_b2n * (*leverarm);
    Z.block<3, 1>(3, 0) = imu_vel - gnss_vel - C_b2n * (*leverarm).cross(gyro) - SkewMat(w_ni_n) * C_b2n * (*leverarm);

	  H.setZero();
    H.block<3, 3>(0, 0) = SkewMat(C_b2n * (*leverarm)); //姿态
	  H.block<3, 3>(0, 3) = I33;  //位置
    H.block<3, 3>(3, 0) = -SkewMat(w_ni_n) * SkewMat(C_b2n * (*leverarm)) - C_b2n * SkewMat((*leverarm).cross(gyro)); //姿态
	  H.block<3, 3>(3, 6) = I33; //速度
    H.block<3, 3>(3, 9) = -SkewMat(C_b2n * (*leverarm)); //陀螺零偏

    // R矩阵
	  R.setZero();
	  R.block<3, 3>(0, 0) = gnss_lla_var.asDiagonal();
	  R.block<3, 3>(3, 3) = gnss_vel_var.asDiagonal();
    //cout << init_P << endl;

	  // 计算出 Kalman 增益矩阵 K
	  MD(23,6) K = init_P * H.transpose() * ((H * init_P * H.transpose() + R).inverse());
	  VD(6) V = Z - H * x;

    x = x + K * V;
	  init_P = (IFF - K * H) * init_P * ((IFF - K * H).transpose()) + K * R * K.transpose();

    V3D Rot_error = x.block<3, 1>(0, 0);
    V3D Pos_error = x.block<3, 1>(3, 0);
    V3D Vel_error = x.block<3, 1>(6, 0);
	  V3D Bias_g_error = x.block<3, 1>(9, 0);
	  V3D Bias_a_error = x.block<3, 1>(12, 0);
    // 位置速度的改正
	  imu_pos(0) = imu_pos(0) - Pos_error(0) * DR_inv(0, 0);
	  imu_pos(1) = imu_pos(1) - Pos_error(1) * DR_inv(1, 1);
	  imu_pos(2) = imu_pos(2) - Pos_error(2) * DR_inv(2, 2);
    imu_vel = gnss_vel - Vel_error;
    init_state.pos = imu_pos;
    init_state.vel = imu_vel;
    // 修改状态变量和方差
    kf_state.change_x(init_state);
    kf_state.change_P(init_P);
  }
}
