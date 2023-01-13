#include "ukf.h"
#include "Eigen/Dense"
#include<iostream>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;



/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_<<1,0,0,0,0,
          0,1,0,0,0,
          0,0,1,0,0,
          0,0,0,1,0,
          0,0,0,0,1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;//1.5//30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;//30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  // State dimension
  n_x_=5;

  // Augmented state dimension
  n_aug_=7;

  // Sigma point spreading parameter
  lambda_=3-n_x_;

  // predicted sigma points matrix
  Xsig_pred_=MatrixXd(n_x_,2*n_aug_+1);
  Xsig_pred_.fill(0.0);
  // Weights of sigma points
  weights_=VectorXd(2*n_aug_+1);
  double weight_0= lambda_/(lambda_+n_aug_);
  double weight = 0.5/((lambda_+n_aug_));
  weights_(0)=weight_0;
  for (int i=0;i<2*n_aug_+1;i++)
  {
    weights_(i)=weight;
  } 
  // Noise matrices
  R_radar_ = MatrixXd(3,3);
  R_radar_.fill(0.0);
  R_radar_(0,0)=std_radr_*std_radr_;
  R_radar_(1,1)=std_radphi_*std_radphi_;
  R_radar_(2,2)=std_radrd_*std_radrd_;

  R_lidar_ = MatrixXd(2,2);
  R_lidar_.fill(0.0);
  R_lidar_(0,0)=std_laspx_*std_laspx_;
  R_lidar_(1,1)=std_laspy_*std_laspy_;

  is_initialized_ = false;
  time_us_ = 0.0;

  NIS_radar_=0;
  NIS_lidar_=0;

}

 UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_)
  {
    if(meas_package.sensor_type_==MeasurementPackage::RADAR)
    {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_d = meas_package.raw_measurements_[2];

      double px = rho*cos(phi);
      double py = rho*sin(phi);
      double vx = rho_d*cos(phi);
      double vy = rho_d*sin(phi);
      
      double v = sqrt(vx*vx + vy*vy);
      x_<<px,py,v,0,0;

    }
    else
    {
      x_<<meas_package.raw_measurements_[0],
          meas_package.raw_measurements_[1],
          0,
          0,
          0;       ;
      // P_<<std_radr_*std_radr_,0,0,0,0,
      //     0,std_radr_*std_radr_,0,0,0,
      //     0,0,std_radrd_*std_radrd_,0,0,
      //     0,0,0,0,0,
      //     0,0,0,0,0; 
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_=true;
    std::cout<<"Inizialization finished"<<std::endl;
    return;
  }
  double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;                                                      
  
  // Prediction step
  Prediction(dt);
  // Update
  if((meas_package.sensor_type_ == MeasurementPackage::LASER))// && use_laser_)
  {
    UpdateLidar(meas_package);
  }
  if((meas_package.sensor_type_ == MeasurementPackage::RADAR))// && use_radar_)
  {
    UpdateRadar(meas_package);
  }
  

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  // Create augmentation
  VectorXd X_aug = VectorXd(n_aug_);
  X_aug.fill(0.0);
  X_aug.head(n_x_)=x_;

  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5)=std_a_*std_a_;
  P_aug(6,6)=std_yawdd_*std_yawdd_;
  //P_aug.bottomRightCorner(2,2)=Q;
  //std::cout<<"P_Aug: "<<P_aug<<std::endl;
  // Create square matrix
  MatrixXd A = P_aug.llt().matrixL();
  // create sigma Point matrix
  MatrixXd Xsig_aug=MatrixXd(n_aug_,2*n_aug_+1);
  Xsig_aug.col(0)=X_aug;
  double factor = sqrt(lambda_+n_aug_);
  for (int column =0;column<n_aug_;column++)
  {
    Xsig_aug.col(column+1) = X_aug+factor*A.col(column);
    Xsig_aug.col(column+1+n_aug_) = X_aug-factor*A.col(column);
  }
  // //std::cout<<"Xsig_aug "<<Xsig_aug<<std::endl;
  // Predict sigma points
  Xsig_pred_.fill(0.0);
  for (int i =0; i< 2 * n_aug_ + 1 ;++i)
  {
    VectorXd sigmaPoint=Xsig_aug.col(i);
    double px       = sigmaPoint(0);
    double py       = sigmaPoint(1);
    double v        = sigmaPoint(2);
    double yaw      = sigmaPoint(3);
    double yawd     = sigmaPoint(4);
    double nu_a     = sigmaPoint(5);
    double nu_yawdd = sigmaPoint(6);

    if (fabs(yawd)>0.001)
    {
      Xsig_pred_(0,i) = px+(v/yawd)*(sin(yaw+delta_t*yawd) - sin(yaw)) + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
      Xsig_pred_(1,i) = py+(v/yawd)*(cos(yaw) - cos(yaw+delta_t*yawd)) + 0.5*delta_t*delta_t*sin(yaw)*nu_a;       
    
    }
    else
    {
      Xsig_pred_(0,i) = px + v*cos(yaw)*delta_t +0.5*delta_t*delta_t*cos(yaw)*nu_a;
      Xsig_pred_(1,i) = py + v*sin(yaw)*delta_t +0.5*delta_t*delta_t*sin(yaw)*nu_a;        
    }
    Xsig_pred_(2,i)=v + delta_t*nu_a;
    Xsig_pred_(3,i)=yaw + yawd*delta_t + 0.5*delta_t*delta_t*nu_yawdd;
    Xsig_pred_(4,i)=yawd + delta_t*nu_yawdd;
        
  }  
  // // Predict state
  // // predicted state mean
  // std::cout<< "X before prediction: "<<x_<<std::endl;
  // VectorXd x_pred=VectorXd(n_x_);
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) >  M_PI) x_diff(3) -=2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) +=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
  std::cout<< "Predicted X: "<<x_<<std::endl; 
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 2;
  VectorXd z = meas_package.raw_measurements_;
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  
  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    double px=Xsig_pred_(0,i);
    double py=Xsig_pred_(1,i);;
    Zsig(0,i)=px;
    Zsig(1,i)=py;
        
  }
  
  // calculate mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  // calculate innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  S=S+ R_lidar_;
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  // calculate cross correlation matrix
  Tc.fill(0.0);
  //std::cout<< "Here! "<<std::endl;
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
     while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
     while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  // calculate Kalman gain K;
  MatrixXd K(n_z,n_z);
  K = Tc*S.inverse();  
  // update state mean and covariance matrix
  VectorXd z_diff=z - z_pred;
  x_ = x_ + K*(z_diff);
  P_ = P_ - K*S*K.transpose();
  // compute normalized innovation squared(NIS)
  NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;
  std::cout<< "Predicted X Lidar: "<<x_<<std::endl;
  std::cout<< "NIS Lidar: "<<NIS_lidar_<<std::endl;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  VectorXd z = meas_package.raw_measurements_;
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // transform sigma points into measurement space
  Zsig.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    double px=Xsig_pred_(0,i);
    double py=Xsig_pred_(1,i);
    double v=Xsig_pred_(2,i);
    double psi=Xsig_pred_(3,i);
    Zsig(0,i)=sqrt(px*px+py*py);
    Zsig(1,i)=atan2(py,px);
    Zsig(2,i)=((px*cos(psi)*v)+(py*sin(psi)*v))/(sqrt(px*px+py*py));    
  }
  // calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  // calculate innovation covariance matrix S  
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
     while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
     while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  S=S+ R_radar_;
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
     while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
     while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

     while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
     while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  // calculate Kalman gain K;
  MatrixXd K(n_z,n_z);
  K = Tc*S.inverse();  
  // update state mean and covariance matrix
  VectorXd z_diff=z - z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
  x_ = x_ + K*(z_diff);
  P_ = P_ - K*S*K.transpose();
  // compute normalized innovation squared(NIS)
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  
  std::cout<< "Predicted X Radar: "<<x_<<std::endl;
  std::cout<< "NIS Radar: "<<NIS_radar_<<std::endl;
}
