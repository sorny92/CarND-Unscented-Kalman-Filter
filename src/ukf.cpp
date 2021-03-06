#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);

  //measurement covariance matrix - laser
  R_laser_ << std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;

  //measurement covariance matrix - radar
  R_radar_ << std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;
}

UKF::~UKF() {}
 
/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: " << endl;
    
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1),0,0,0;
    } else {
      const double rho = meas_package.raw_measurements_(0);
      const double phi = meas_package.raw_measurements_(1);
      const double d_rho = meas_package.raw_measurements_(2);

      const double module_d_rho = sqrt(d_rho*cos(phi)*d_rho*cos(phi) +
                                       d_rho*sin(phi)*d_rho*sin(phi));
      x_ << rho * cos(phi),
            rho * sin(phi),
            module_d_rho,
            0,
            0;
    }
    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ =meas_package.timestamp_;
    return;
  }

  //Convert to seconds
  double delta_t = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;
  cout<< "delta_t: " << delta_t << endl;
  /*while (delta_t > 0.1)
  {
  
  //  cout << "Press ENTER to continue...";
  //  cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
  
  const double dt = 0.05;
  Prediction(dt);
  delta_t -= dt;
  }*/

          
  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR ) {
    // Radar updates
    cout << " RADAR" << endl;
    UpdateRadar(meas_package);
  } else {
    // Laser updates
    cout << " LIDAR" << endl;
    UpdateLidar(meas_package);
  }

  previous_timestamp_ = meas_package.timestamp_;
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
  /*************************************************************************
   * GENERATION SIGMA POINTS
  **************************************************************************/
  
  MatrixXd X_sig_aug_ = MatrixXd(n_aug_, n_sigma_points);

  VectorXd x_aug_ = VectorXd(n_aug_);
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);

  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;

  MatrixXd A = P_aug_.llt().matrixL();

  // Column 0 equal to the mean
  X_sig_aug_.col(0) = x_aug_;
  // Remaining columns following the formula
  for (int i = 0; i< n_aug_; i++){
    X_sig_aug_.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_) * A.col(i);
    X_sig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  /*************************************************************************
   * PREDICT SIGMA POINTS
  *************************************************************************/
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_points);
  for (int i = 0; i< n_sigma_points; i++)
  {
    //extract values for better readability
    double p_x = X_sig_aug_(0,i);
    double p_y = X_sig_aug_(1,i);
    double v = X_sig_aug_(2,i);
    double yaw = X_sig_aug_(3,i);
    double yawd = X_sig_aug_(4,i);
    double nu_a = X_sig_aug_(5,i);
    double nu_yawdd = X_sig_aug_(6,i);

    //predicted state values
    double px_p =0, py_p=0;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
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
  std::cout << "X_sig: " << std::endl << Xsig_pred_ << std::endl;
  /*******************************************************************
   * USE SIGMA POINTS TO CALCULATE MEAN AND COVARIANCE
  *******************************************************************/
  weights_ = VectorXd::Zero(n_sigma_points);
  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i< n_sigma_points; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sigma_points; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sigma_points; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = atan2(sin(x_diff(3)),cos(x_diff(3)));
    //or
    /*while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;*/

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
  //mean predicted measurement
  z_ = VectorXd(2);
  z_[0] = meas_package.raw_measurements_[0];
  z_[1] = meas_package.raw_measurements_[1];

  
  MatrixXd z_sigma_points = MatrixXd(2, n_sigma_points);
  for (int i = 0; i < n_sigma_points; i++) {  //2n+1 sigma points
        // measurement model
        z_sigma_points(0,i) = Xsig_pred_(0,i);          //px
        z_sigma_points(1,i) = Xsig_pred_(1,i);          //py
        
  }
  
  //mean predicted measurement
  z_pred = VectorXd(2);
  z_pred.fill(0.0);
  for (int i=0; i < n_sigma_points; i++) {
      z_pred = z_pred + weights_(i) * z_sigma_points.col(i);
  }


  S_ = MatrixXd(2,2);
  S_.fill(0.0);

  MatrixXd Tc_ = MatrixXd(5,2);
  Tc_.fill(0.0);


  for (int i = 0; i < n_sigma_points; i++) {  //2n+1 simga points
    VectorXd z_diff = VectorXd(2);
    z_diff = z_sigma_points.col(i) - z_pred;
    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc_ = Tc_ + weights_(i) * x_diff * z_diff.transpose();
  }
  //add measurement noise covariance matrix
  S_ = S_ + R_laser_;

  /**************************************************************************
 * UPDATE STEP
 *************************************************************************/

  MatrixXd K = Tc_ * S_.inverse();

  VectorXd z_diff = VectorXd::Zero(2);
  z_diff = z_ - z_pred;

    //update state mean and covariance matrix
  //cout << K << "\n" << endl;
  //cout << z_diff << "\n" << endl;
  x_ = x_ + K * z_diff;
  //cout << x_ << "\n" << endl;
  P_ = P_ - K*S_*K.transpose();
  //cout << P_ << "\n\n" << endl;

  NIS_laser_ = z_diff.transpose() * S_.inverse() * z_diff;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  z_ = VectorXd(3);
  z_[0] = meas_package.raw_measurements_[0];
  z_[1] = meas_package.raw_measurements_[1];
  z_[2] = meas_package.raw_measurements_[2];
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //
  MatrixXd z_sigma_points = MatrixXd(3, n_sigma_points);

  //transform sigma points into measurement space
  for (int i = 0; i < n_sigma_points; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);


    if(fabs(p_x) <= 0.001){
      p_x = 0.001;
    }
    if(fabs(p_y) <= 0.001){
      p_y = 0.001;
    }
    // measurement model
    const double r = sqrt(p_x*p_x + p_y*p_y);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    z_sigma_points(0,i) = r;                        //r
    z_sigma_points(1,i) = atan2(p_y,p_x);           //phi 
    z_sigma_points(2,i) = (p_x*v1 + p_y*v2 ) / r;   //r_dot
  }
  //cout << z_sigma_points << "\n" << endl;
  //mean predicted measurement
  z_pred = VectorXd(3);
  z_pred.fill(0.0);
  for (int i=0; i < n_sigma_points; i++) {
      z_pred = z_pred + weights_(i) * z_sigma_points.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S_ = MatrixXd(3,3);
  S_.fill(0.0);

  MatrixXd Tc_ = MatrixXd(5,3);
  Tc_.fill(0.0);

  
  for (int i = 0; i < n_sigma_points; i++) {  //2n+1 simga points
    VectorXd z_diff = VectorXd(3);
    //residual
     z_diff = z_sigma_points.col(i) - z_pred;

    //angle normalization
    if (z_diff(1)> M_PI){
      z_diff(1) = fmod(z_diff(1),(2.*M_PI)) - 2.*M_PI;
    }else if(z_diff(1)< -M_PI){
      z_diff(1) = fmod(z_diff(1),(2.*M_PI)) + 2.*M_PI;
    }

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    if (x_diff(3)> M_PI){
      x_diff(3) = fmod(x_diff(3),(2.*M_PI)) - 2.*M_PI;
    }else if(x_diff(3)< -M_PI){
      x_diff(3) = fmod(x_diff(3),(2.*M_PI)) + 2.*M_PI;
    }

    Tc_ = Tc_ + weights_(i) * x_diff * z_diff.transpose();
  }
  //cout << z_diff << "\n" << endl;
  //add measurement noise covariance matrix
  S_ = S_ + R_radar_;

/**************************************************************************
 * UPDATE STEP
 *************************************************************************/
  //Kalman gain K;
  // TODO: Tc_ and S_ have -nan or really close to 0 values (1e-300)
  //cout << Tc_ << "\n" << endl;
  //cout << S_ << "\n" << endl;
  MatrixXd K = Tc_ * S_.inverse();

  //residual
  //TODO: Revisar z_pred
  //cout << z_pred << "\n" << endl;
  VectorXd z_diff = VectorXd(3);
  z_diff = z_ - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  //cout << K << "\n" << endl;
  //cout << z_diff << "\n" << endl;
  x_ = x_ + K * z_diff;
  //cout << x_ << "\n" << endl;
  P_ = P_ - K*S_*K.transpose();
  //cout << P_ << "\n\n" << endl;

  NIS_radar_ = z_diff.transpose() * S_.inverse() * z_diff;
}
