#include <Eigen/Eigen>

#include "use-ikfom.hpp"
#include <typeinfo>
#include <GeographicLib/LocalCartesian.hpp>

using namespace GeographicLib;
/*============================================定义一些地理坐标常数==============================================*/
/* Physical parameters of the Earth, Sun and Moon  */
#define R_WGS84 6378137.0           /* Radius Earth [m]; WGS-84  */
#define E2_WGS84 0.0066943799901413
#define Omega_WGS  7.2921151467e-5    /*[rad/s], the earth rotation rate */

/*===========================================================================================================*/

/*******************************************************
 * @brief       计算卯酉圈曲率半径
 * @author      Sihan Wang
 * @paras[in]   fai（纬度，rad）
 * @paras[out]  none
 * @return      RN(double)
 * @Date        2022/5/20
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
 * @return      RM(double)
 * @Date        2022/5/20
 *******************************************************/
double cal_RM(double fai) {
  double W = sqrt(1.0 - E2_WGS84 * sin(fai) * sin(fai));
  double RM = R_WGS84 * (1.0 - E2_WGS84) / (W * W * W);
  return RM;
}

/*******************************************************
 * @brief       反对称矩阵
 * @author      Sihan Wang
 * @paras[in]   Vector3d
 * @paras[out]  none
 * @return      反对称矩阵Skew(Matrix3d)
 * @Date        2022/5/20
 *******************************************************/
M3D SkewMat(V3D Vec) {
  M3D Skew = Matrix3d::Zero(3, 3);
  Skew(0, 1) = -Vec(2);
  Skew(0, 2) = Vec(1);
  Skew(1, 0) = Vec(2);
  Skew(1, 2) = -Vec(0);
  Skew(2, 0) = -Vec(1);
  Skew(2, 1) = Vec(0);
  return Skew;
}

/*******************************************************
 * @brief       旋转矩阵标准化
 * @author      Kaiserreich-Unofficial
 * @paras[in]   Matrix3d
 * @paras[out]  none
 * @return      Matrix3d
 * @Date        2024/3/31
 *******************************************************/
M3D NormalizeR(const M3D &R)
{
    JacobiSVD<M3D> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);
    M3D U = svd.matrixU();
    M3D V = svd.matrixV();
    M3D S = Matrix3d::Identity();
    S(2, 2) = (U * V.transpose()).determinant();
    return U * S * V.transpose();
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
    void CorrectLeverArm(const V3D &leverarm, const V3D &gyro);

    /* 姿态变量 */
    bool gnss_need_init_ = true;
    V3D imu_pos;
    V3D imu_vel;
    /* 定义地理坐标变量 */
    // 创建LocalCartesian对象，用于ENU和ECEF之间的转换
    LocalCartesian proj;

    double lat, lon, alt;
    double vE, vN, vU;
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
    : gnss_need_init_(true)
{}

GNSSProcess::~GNSSProcess() {}

void GNSSProcess::GNSS_Init(const V3D &init_pos) {
  if(init_pos.sum() <= 1e-3)
  {
    ROS_WARN("GPS NavSat Unfixed!");
  }
  else
  {
    proj.Reset(init_pos(0), init_pos(1), init_pos(2));
    gnss_need_init_ = false;
    ROS_INFO("GNSS Initial Done");
  }
}
void GNSSProcess::CorrectLeverArm(const V3D &leverarm, const V3D &gyro)
{
  M3D mat_n2l = Matrix3d::Zero(3, 3); //将北东地转化为LLH系

	mat_n2l(0, 0) = 1.0 / ((RN + alt) * cos(lat));
  mat_n2l(1, 1) = 1.0 / (RM + alt);
	mat_n2l(2, 2) = 1.0;

  M3D OMIGA_ni_n = SkewMat(w_ni_n);
	M3D OMIGA_bi_b = SkewMat(gyro);
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

  V3D imu_enu_pos = init_state.pos;
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
    ROS_INFO("[IMU]: E:%0.6f N:%0.6f U:%0.6f", imu_enu_pos(0), imu_enu_pos(1), imu_enu_pos(2));
    // 将ENU坐标转换为LLA坐标
    proj.Reverse(imu_enu_pos(0), imu_enu_pos(1), imu_enu_pos(2), imu_pos(0), imu_pos(1), imu_pos(2));
    ROS_INFO("[IMU]: Lat:%0.6f Lon:%0.6f Alt:%0.6f", imu_pos(0), imu_pos(1), imu_pos(2));

    lat = imu_pos(0);
    lon = imu_pos(1);
    alt = imu_pos(2);

    vE = imu_vel(0);
    vN = imu_vel(1);
    vU = imu_vel(2);

    RM = cal_RM(lat);
    RN = cal_RN(lat);

    w_ei_n << 0, Omega_WGS * cos(lat), Omega_WGS * sin(lat);
    w_ne_n << -1.0 * vN / (RM + alt), vE / (RN + alt), vE * tan(lat) / (RN + alt);
    w_ni_n = w_ei_n + w_ne_n;

    C_b2n = init_state.rot.toRotationMatrix();

    ROS_INFO("[GNSS]: Lat:%0.6f Lon:%0.6f Alt:%0.6f", gnss_lla_cord(0), gnss_lla_cord(1), gnss_lla_cord(2));

    // CorrectLeverArm(*leverarm, gyro); //杆臂效应改正

    VD(23) x;
    x.setZero();
    // x.block<3, 1>(0, 0) = SO3ToEuler(init_state.rot);
    x.block<3, 1>(0, 3) = imu_pos;
    x.block<3, 1>(0, 6) = init_state.vel;
    // x.block<3, 1>(0, 9) = init_state.bg;
    // x.block<3, 1>(0, 12) = init_state.ba;

    // DR_inv 将 n 系下的北向、东向和垂向位置差异（单位 m）转化为纬度、经度和高程分量的差异
    M3D DR;
    DR.setZero();
    DR(0, 0) = (RN + alt) * cos(lat);
    DR(1, 1) = RM + alt;
    DR(2, 2) = 1;
    M3D DR_inv = DR.inverse();

    /***********************观测方程***************************/
    Z.setZero();
    Z.block<3, 1>(0, 0) = DR * (imu_pos - gnss_pos) + C_b2n * (*leverarm);
    Z.block<3, 1>(3, 0) = imu_vel - gnss_vel - C_b2n * (*leverarm).cross(gyro) - SkewMat(w_ni_n) * C_b2n * (*leverarm);

	  H.setZero();
    // H.block<3, 3>(0, 0) = SkewMat(C_b2n * (*leverarm)); //姿态
	  H.block<3, 3>(0, 3) = I33;  //位置
    // H.block<3, 3>(3, 0) = -SkewMat(w_ni_n) * SkewMat(C_b2n * (*leverarm)) - C_b2n * SkewMat((*leverarm).cross(gyro)); //姿态
	  H.block<3, 3>(3, 6) = I33; //速度
    // H.block<3, 3>(3, 9) = -SkewMat(C_b2n * (*leverarm)); //陀螺零偏

    // R矩阵
	  R.setZero();
	  R.block<3, 3>(0, 0) = gnss_lla_var.asDiagonal();
	  R.block<3, 3>(3, 3) = gnss_vel_var.asDiagonal();
    // cout << init_P << endl;

	  // 计算出 Kalman 增益矩阵 K
	  MD(23,6) K = init_P * H.transpose() * ((H * init_P * H.transpose() + R).inverse());
	  VD(6) V = Z - H * x;

    x = x + K * V;
	  init_P = (IFF - K * H) * init_P * ((IFF - K * H).transpose()) + K * R * K.transpose();

    // V3D _rot_error = x.block<3, 1>(0, 0);
    V3D _pos_error = x.block<3, 1>(3, 0);
    V3D _vel_error = x.block<3, 1>(6, 0);
	  // V3D _bias_g_error = x.block<3, 1>(9, 0);
	  // V3D _bias_a_error = x.block<3, 1>(12, 0);
    // 位置速度的改正
	  imu_pos(0) -= _pos_error(0) * DR_inv(0, 0);
	  imu_pos(1) -= _pos_error(1) * DR_inv(1, 1);
	  imu_pos(2) -= _pos_error(2) * DR_inv(2, 2);
    imu_vel = gnss_vel - _vel_error;
    /* 将改正后的位置和速度写回init_state中 */
    // 将NED坐标转换为LLA坐标
    V3D imu_temp;
    proj.Forward(imu_pos(0), imu_pos(1), imu_pos(2), imu_temp(0), imu_temp(1), imu_temp(2));
    imu_pos = imu_temp;

    init_state.pos = imu_pos;
    init_state.vel = imu_vel;
    // 姿态改正
    // M3D _C_b2n_new = (I33 - SkewMat(_rot_error)).inverse() * C_b2n;
    // init_state.rot = NormalizeR(_C_b2n_new);

    // 角加速度零偏改正
    // init_state.bg += _bias_g_error;
    // init_state.ba += _bias_a_error;
    // 修改状态变量和方差
    kf_state.change_x(init_state);
    kf_state.change_P(init_P);
  }
}
