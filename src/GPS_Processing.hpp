#include <Eigen/Eigen>

#include "use-ikfom.hpp"

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

    GNSSProcess();
    ~GNSSProcess();
    Eigen::Matrix<double, 6, 1> Z;  // Z: 观测矩阵(pos+vel)
    deque<GNSSPtr> v_gnss;

  void Process(const MeasureGroup &meas,
               esekfom::esekf<state_ikfom, 12, input_ikfom> &kf_state);

  private:
    GNSSPtr last_gnss_;
    bool init = false;
    V3D init_pos;

    void GNSS_Init(const V3D &lla_cord);
    V3D convertToEnu(const V3D &lla_cord);
    void Update(const Eigen::VectorXd &v_gnss);
};

GNSSProcess::GNSSProcess() { GNSSPtr last_gnss_ = new Eigen::VectorXd(13); }

GNSSProcess::~GNSSProcess() {}

V3D GNSSProcess::convertToEnu(const V3D &lla_cord) {
  V3D enu;
  enu.setZero();
  if (!init) {
    ROS_WARN("GNSS Convert To ENU Cord Function Uninited!");
    return enu;
  }

  double radLat1 = deg2rad(init_pos(0));
  double radLong1 = deg2rad(init_pos(1));
  double radLat2 = deg2rad(lla_cord(0));
  double radLong2 = deg2rad(lla_cord(1));

  double delta_lat = radLat2 - radLat1;
  double delta_long = radLong2 - radLong1;

  enu(1) = R_WGS84 * delta_lat;
  enu(0) = R_WGS84 * cos(radLat1) * delta_long;
  enu(2) = lla_cord(2) - init_pos(2);
  ROS_INFO("x:%0.6f y:%0.6f z:%0.6f", enu(0), enu(1), enu(2));

  return enu;
}

void GNSSProcess::GNSS_Init(const V3D &lla_cord) {
  init_pos = lla_cord;
  init = true;
  ROS_INFO("GNSS Initial Done");
}

void GNSSProcess::Update(const Eigen::VectorXd &v_gnss) {
  const double &gnss_time = v_gnss(0);
  const V3D &gnss_neu_cord = convertToEnu(v_gnss.segment(1, 3));
  const V3D &gnss_vel = v_gnss.segment(7, 9);
}

void GNSSProcess::Process(const MeasureGroup &meas, esekfom::esekf<state_ikfom, 12, input_ikfom> &kf_state) {
  /*** add the gnss of the last frame-tail to the of current frame-head ***/
  deque<GNSSPtr> v_gnss = meas.gnss;
  // v_gnss.push_front(last_gnss_);
  /* 获取当前系统状态变量 */
  state_ikfom imu_state = kf_state.get_x();
  const V3D &imu_pos = imu_state.pos;
  const V3D &imu_vel = imu_state.vel;
  cout << imu_pos << imu_vel << endl;
  // V3D imu_eul = SO3ToEuler(imu_state.rot);
  /* 获取当前系统状态变量 */
  if (!init) {
    const V3D &gnss_lla_cord = v_gnss.front()->segment(1, 3);
    cout << gnss_lla_cord << endl;
    GNSS_Init(gnss_lla_cord);
  }

  auto gnss_iter = v_gnss.begin();
  while (gnss_iter != v_gnss.end()) {
    Update(*(v_gnss.front()));
    gnss_iter = v_gnss.erase(gnss_iter);
  }
  cout << "GNSS DATA Processed Done!";
  // cout << imu_pose << endl;
  // cout << imu_vel << endl;
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

  M3D C_b2n = imu_state.rot.toRotationMatrix();

  // DR_inv 将 n 系下的北向、东向和垂向位置差异（单位
  // m）转化为纬度、经度和高程分量的差异
  M3D DR;
  DR.setZero();
  DR(0, 0) = RM + h;
  DR(1, 1) = (RN + h) * cos(B);
  DR(2, 2) = -1;
  M3D DR_inv = DR.inverse();

  /***********************观测方程***************************/
  Z.setZero();
  // Z.block<3, 1>(0, 0) = DR * (imu_pos - gnss_pos);
  // Z.block<3, 1>(3, 0) = imu_vel - gnss_vel;
}
