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

    void Process(const MeasureGroup &meas, esekfom::esekf<state_ikfom, 12, input_ikfom> &kf_state);

    V3D mean_lla_cord;
    V3D mean_lla_var;
    V3D mean_vel;
    V3D mean_vel_var;

  private:
    void GNSS_Init(const V3D &lla_cord);
    V3D ConvertToNED(const V3D &lla_cord);

    bool gnss_need_init_ = true;
    V3D* init_pos;
    M3D I33 = Matrix3d::Identity();
    MD(23,23) IFF = MatrixXd::Identity(23,23);
    // void Update(const Eigen::VectorXd &v_gnss);
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

void GNSSProcess::Process(const MeasureGroup &meas, esekfom::esekf<state_ikfom, 12, input_ikfom> &kf_state) {
  /* 获取当前系统状态变量 */
  state_ikfom init_state = kf_state.get_x(); //rot pos vel bias_gyro bias_acc gravity R_L p_L
  /* 获取当前系统协方差矩阵 */
  esekfom::esekf<state_ikfom, 12, input_ikfom>::cov init_P = kf_state.get_P(); // 24x24维

  V3D imu_pos = init_state.pos;
  V3D imu_vel = init_state.vel;
  // V3D imu_eul = SO3ToEuler(init_state.rot);
  /* 获取当前系统状态变量 */
  uint8_t N = 1;
  for (const auto &gnss : meas.gnss)
  {
    const auto &lla_cord = gnss->segment<3>(1);
    const auto &lla_var = gnss->segment<3>(4);
    const auto &vel = gnss->segment<3>(7);
    const auto &vel_var = gnss->segment<3>(10);

    mean_lla_cord += (lla_cord - mean_lla_cord) / N;
    mean_lla_var  += (lla_var - mean_lla_var) / N;
    mean_vel      += (vel - mean_vel) / N;
    mean_vel_var  += (vel_var - mean_vel_var) / N;
    // cout<<"acc norm: "<<cur_acc.norm()<<" "<<mean_acc.norm()<<endl;
    N++;
  }
  /* 若GNSS需要初始化 */
  if(gnss_need_init_) GNSS_Init(mean_lla_cord);
  else
  {
    double B = imu_pos(0);
    double L = imu_pos(1);
    double h = imu_pos(2);

    double vN = imu_vel(0);
    double vE = imu_vel(1);
    double vD = imu_vel(2);

    double RM = cal_RM(B);
    double RN = cal_RN(B);

    V3D w_ei_n;  // 地球自转角速度在n系下的投影
    V3D w_ne_n;  // n系相对于e系的角速度在n系下的投影（牵连角速度）
    V3D w_ni_n;  // n系相对于i系的角速度在n系下的投影（牵连角速度）
    w_ei_n << Omega_WGS * cos(B), 0, -1.0 * Omega_WGS * sin(B);
    w_ne_n << vE / (RN + h), -1.0 * vN / (RM + h), -1.0 * vE * tan(B) / (RN + h);
    w_ni_n = w_ei_n + w_ne_n;

    M3D C_b2n = init_state.rot.toRotationMatrix();
    VD(23) x;
    x.setZero();
    x.block<3, 1>(0, 3) = init_state.pos;
    x.block<3, 1>(0, 6) = init_state.vel;

    // DR_inv 将 n 系下的北向、东向和垂向位置差异（单位 m）转化为纬度、经度和高程分量的差异
    M3D DR;
    DR.setZero();
    DR(0, 0) = RM + h;
    DR(1, 1) = (RN + h) * cos(B);
    DR(2, 2) = -1;
    M3D DR_inv = DR.inverse();

    /***********************观测方程***************************/
    V3D gnss_pos = ConvertToNED(mean_lla_cord);
    ROS_INFO("Lat:%0.6f Lon:%0.6f Alt:%0.6f", mean_lla_cord(0), mean_lla_cord(1), mean_lla_cord(2));
    ROS_INFO("N:%0.6f E:%0.6f D:%0.6f", gnss_pos(0), gnss_pos(1), gnss_pos(2));
    Z.setZero();
    Z.block<3, 1>(0, 0) = DR * (imu_pos - gnss_pos); // 还没有加入leverarm效应修正
    Z.block<3, 1>(3, 0) = imu_vel - mean_vel;

	  H.setZero();
	  H.block<3, 3>(0, 3) = I33;
	  H.block<3, 3>(3, 6) = I33;

    // R矩阵
	  R.setZero();
	  R.block<3, 3>(0, 0) = mean_lla_var.asDiagonal();
	  R.block<3, 3>(3, 3) = mean_vel_var.asDiagonal();
    //cout << init_P << endl;

	  // 计算出 Kalman 增益矩阵 K
	  MD(23,6) K = init_P * H.transpose() * ((H * init_P * H.transpose() + R).inverse());
	  VD(6) V = Z - H * x;

    x = x + K * V;
	  init_P = (IFF - K * H) * init_P * ((IFF - K * H).transpose()) + K * R * K.transpose();

    V3D Pos_error = x.block<3, 1>(0, 3);
    V3D Vel_error = x.block<3, 1>(0, 6);
    // 位置速度的改正
	  imu_pos(0) = imu_pos(0) - Pos_error(0) * DR_inv(0, 0);
	  imu_pos(1) = imu_pos(1) - Pos_error(1) * DR_inv(1, 1);
	  imu_pos(2) = imu_pos(2) - Pos_error(2) * DR_inv(2, 2);
    imu_vel = mean_vel - Vel_error;
    init_state.pos = imu_pos;
    init_state.vel = imu_vel;
    kf_state.change_x(init_state);
  }
}
