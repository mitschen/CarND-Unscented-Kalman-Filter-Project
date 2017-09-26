#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
  : is_initialized_(false)
  , use_laser_(true)
  , use_radar_(true)
  , n_x_(5)   //[pos1 pos2 vel_abs yaw_angle yaw_rate]
  , n_aug_(7) //[pos1 pos2 vel_abs yaw_angle yaw_rate, noise_vel, noise_yaw]
  , lambda_(3-n_aug_)
  , n_sigma_pts_(2*n_aug_+1)
  , x_(n_x_)
  , P_(n_x_,n_x_)
  , Xsig_pred_(n_x_, n_sigma_pts_)
  , Q_(2,2)
  , R_radar_(3,3)
  , R_laser_(2,2)
  , time_us_(0LL)
  , std_a_(30)
  , std_yawdd_(30)
  , std_laspx_(0.15)
  , std_laspy_(0.15)
  , std_radr_(0.3)
  , std_radphi_(0.03)
  , std_radrd_(0.3)
  , weights_(n_sigma_pts_)
  {



  // initial state vector
  // will be initialized with the first measurement values
  x_.setZero();

//  // initial covariance matrix
//  P_.setIdentity();
  P_ <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  //set the sigma ptr prediction to zero
  Xsig_pred_.setZero();

//  // Process noise standard deviation longitudinal acceleration in m/s^2
//  std_a_ = 30;
//  // Process noise standard deviation yaw acceleration in rad/s^2
//  std_yawdd_ = 30;
//  // Laser measurement noise standard deviation position1 in m
//  std_laspx_ = 0.15;
//  // Laser measurement noise standard deviation position2 in m
//  std_laspy_ = 0.15;
//  // Radar measurement noise standard deviation radius in m
//  std_radr_ = 0.3;
//  // Radar measurement noise standard deviation angle in rad
//  std_radphi_ = 0.03;
//  // Radar measurement noise standard deviation radius change in m/s
//  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  Q_ << std_a_ * std_a_,0.,
        0., std_yawdd_*std_yawdd_;
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0, std_radrd_*std_radrd_;
  R_laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

  //initialize the weights
  weights_.fill(0.5/(n_aug_+lambda_));
  weights_(0) = lambda_/(lambda_+n_aug_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
   double const dt( (meas_package.timestamp_ - time_us_) / 1000000.0);
   if(use_radar_ && MeasurementPackage::RADAR == meas_package.sensor_type_)
   {
      double const &rho(meas_package.raw_measurements_[0]);
      double const &phi(meas_package.raw_measurements_[1]);
      double const &rhodot(meas_package.raw_measurements_[2]);
      double const x(cos(phi) * rho);
      double const y(sin(phi) * rho);
      if(!is_initialized_)
      {
         x_ << x, y, 0.41 /*15km/h*/, 0, 0;
         time_us_=meas_package.timestamp_;
         is_initialized_ = true;
         return;
      }

      Prediction(dt);
      UpdateRadar(meas_package);
   }
   else if(use_laser_)
   {
      double const &x(meas_package.raw_measurements_[0]);
      double const &y(meas_package.raw_measurements_[1]);
      if(!is_initialized_)
      {
         x_ << x, y, 0.41 /*15km/h*/, 0, 0;
         time_us_=meas_package.timestamp_;
         is_initialized_ = true;
         return;
      }
      Prediction(dt);
      UpdateLidar(meas_package);
   }
   time_us_=meas_package.timestamp_;

   std::cout<<"StateVectorX "<<std::endl<<x_<<std::endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */


  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  ///_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_
  //GENERATE THE SIGMA POINTS
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug.bottomRightCorner(2,  2) = Q_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  MatrixXd Xsig_aug(n_aug_, n_sigma_pts_);
  Xsig_aug.col(0)  = x_aug;
  static double const sqrtLamda_nAug(sqrt(lambda_+n_aug_));
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrtLamda_nAug * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrtLamda_nAug * L.col(i);
  }

  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  ///_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_
  //PREDICT THE SIGMA POINTS
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  //now apply sigma point prediction base on the augmented sigma points
  for (int i(0); i< n_sigma_pts_; i++)
  {
    //extract values for better readability
    double const p_x(Xsig_aug(0,i));
    double const p_y(Xsig_aug(1,i));
    double const v(Xsig_aug(2,i));
    double const yaw(Xsig_aug(3,i));
    double const yawd(Xsig_aug(4,i));
    double const nu_a(Xsig_aug(5,i));
    double const nu_yawdd(Xsig_aug(6,i));

    //predicted state values
    double px_p, py_p;

    double const sinyaw(sin(yaw));
    double const cosyaw(cos(yaw));
    double const yawdelta_t(yawd*delta_t);
    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      double const vDyawd(v/yawd);
       px_p = p_x + vDyawd * ( sin (yaw + yawdelta_t) - sinyaw);
       py_p = p_y + vDyawd * ( cosyaw - cos(yaw+yawdelta_t) );
    }
    else {
      double const vdelta_t(v*delta_t);
       px_p = p_x + vdelta_t*cosyaw;
       py_p = p_y + vdelta_t*sinyaw;
    }

    double v_p = v;
    double yaw_p = yaw + yawdelta_t;
    double yawd_p = yawd;

    double const nu_a_delta_t2_05(0.5*nu_a*delta_t*delta_t);

    //add noise
    px_p = px_p + nu_a_delta_t2_05 * cosyaw;
    py_p = py_p + nu_a_delta_t2_05 * sinyaw;
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
  ///_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_
  //PREDICT MEAN AND COVARIANCE
  //_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

  //start with the predicted mean state first
  x_.fill(0.);
  for (int i(0); i < n_sigma_pts_; i++) {
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //and now identify the covariance matrix
  P_.fill(0.);
  for (int i(0); i < n_sigma_pts_; i++) {

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //for lidar we'll only get px, py
    int const n_z(2);

    //current measurement
    VectorXd z = meas_package.raw_measurements_;
    //sanity check
    assert(z.size() == n_z);

    MatrixXd S(n_z, n_aug_);
    S.fill(0.);
    //mean predicted measurement
    VectorXd z_pred(n_z);
    z_pred.fill(0.0);

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, n_sigma_pts_);

    //transform sigma points into measurement space
    //and calculate the mean
    for (int i(0); i < n_sigma_pts_; i++) {

      // extract values for better readibility
      double p_x = Xsig_pred_(0,i);
      double p_y = Xsig_pred_(1,i);

      // measurement model
      Zsig(0,i) = p_x;
      Zsig(1,i) = p_y;

      //calculate the mean prediction measurement
      z_pred = z_pred + weights_(i) * Zsig.col(i);
    }


    //measurement covariance matrix S
    //initialize by the measurement noise covariance matrix
    S = R_laser_;
    for (int i(0); i < n_sigma_pts_; i++) {
      //residual
      VectorXd z_diff = Zsig.col(i) - z_pred;

      S = S + weights_(i) * z_diff * z_diff.transpose();
    }



    //cross correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.);
    for (int i(0); i < n_sigma_pts_; i++) {  //2n+1 simga points

      //residual
      VectorXd z_diff = Zsig.col(i) - z_pred;

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;

      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    MatrixXd const S_inverse = S.inverse();
    //Kalman gain K;
    MatrixXd K = Tc * S_inverse;

    //residual
    VectorXd z_diff = z - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();


    double const NIS(z_diff.transpose()*S_inverse*z_diff);
    std::cout << "NIS (laser): " << NIS << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //for radar we have 3 states
  int const n_z(3);

  //current measurement
  VectorXd z = meas_package.raw_measurements_;
  //sanity check
  assert(z.size() == n_z);

  MatrixXd S(n_z, n_z);
  S.fill(0.);
  //mean predicted measurement
  VectorXd z_pred(n_z);
  z_pred.fill(0.0);

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma_pts_);

  //transform sigma points into measurement space
  //and calculate the mean
  for (int i(0); i < n_sigma_pts_; i++) {

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);          //r
    Zsig(1,i) = atan2(p_y,p_x);                   //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / Zsig(0,i);   //r_dot
    //calculate the mean prediction measurement
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }


  //measurement covariance matrix S
  //initialize by the measurement noise covariance matrix
  S = R_radar_;
  for (int i(0); i < n_sigma_pts_; i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }



  //cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.);
  for (int i(0); i < n_sigma_pts_; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd const S_inverse = S.inverse();
  //Kalman gain K;
  MatrixXd K = Tc * S_inverse;

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();


  double const NIS(z_diff.transpose()*S_inverse*z_diff);

  std::cout << "NIS (radar): " << NIS << std::endl;
}
