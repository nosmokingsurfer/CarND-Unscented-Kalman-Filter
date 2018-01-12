#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

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

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.245;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .625;

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
  TODO:
  
  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  x_.fill(0.0);
  P_.fill(0.0);
  Xsig_pred_ = MatrixXd(5,15);
  Xsig_pred_.fill(0.0);
  time_us_ = 0;

  
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  NIS_radar_ = 0;
  NIS_laser_ = 0;

  weights_ = VectorXd(2*n_aug_ + 1);
  weights_.fill(0.0);
  weights_[0] = lambda_/(lambda_ + n_aug_);
  //cout << "weights_[" << 0 << "] = " << weights_[0] << endl;
  for(int i = 1; i < 2*n_aug_ + 1; i++)
  {
    weights_[i] = 1/(2*(lambda_ + n_aug_));
    //cout << "weights_[" << i << "] = " << weights_[i] << endl;
  }

  //cout << "Weights:" << endl;
  for(int i = 0; i < 2*n_aug_ + 1; i++)
  {
   // cout << "weights_[" << i << "] = " << weights_[i] << endl;
  }

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

  if (!is_initialized_)
  {
    /**
    Initialize the mean x
    */
    double px = 0;
    double py = 0;
    double v = 0;
    double phi = 0;
    double phi_d = 0;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double rho = meas_package.raw_measurements_[0];
      double azim = meas_package.raw_measurements_[1];
      double rho_d = meas_package.raw_measurements_[2];

      px = rho*cos(azim);
      py = rho*sin(azim);
      v = rho_d;
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      px = meas_package.raw_measurements_[0];
      py = meas_package.raw_measurements_[1];
    }

    x_ << px, py, v, phi, phi_d;
    P_ << 0.1,0,0,0,0,
          0,0.1,0,0,0,
          0,0,0.1,0,0,
          0,0,0,0.1,0,
          0,0,0,0,0.1;

    /**
    Initialize the covariance matrix P
    */

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
  }
  else
  {
    double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
    time_us_ = meas_package.timestamp_;

    //cout << delta_t << endl;
    if (delta_t >= 0)
    {
      Prediction(delta_t);
    }

    if((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_)
    {
      UpdateLidar(meas_package);
    }
    else if((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_)
    {
      UpdateRadar(meas_package);
    }
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  //if (delta_t == 0)
  //{
   // cout << delta_t << endl;
    //return;
  //}

  //cout << delta_t << endl;

  MatrixXd Xsig_aug(n_aug_, 2*n_aug_ + 1);
  Xsig_aug.fill(0.0);


  //cout << "Xsig_aug before = " << endl << Xsig_aug << endl;
  GenerateSigmaPoints(&Xsig_aug);
  //cout << "Xsig_aug after = " << endl << Xsig_aug << endl;


  //Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  for(int i = 0; i < 2*n_aug_ + 1; i++)
  {
    VectorXd x_aug = Xsig_aug.col(i);
    VectorXd x = x_aug.head(n_x_);

    double px = x[0];
    double py = x[1];
    double v = x[2];
    double psi = x[3];
    double psi_d = x[4];
    double nu_a = x_aug[5];
    double nu_psidd = x_aug[6];

    double dt_2 = delta_t*delta_t;
    //cout << "dt = " << delta_t << endl;
    //cout << "dt_2 = " << dt_2 << endl;


    VectorXd F(n_x_);
    F.fill(0.0);

    if (fabs(psi_d) > 0.01)
    {
      F << v/psi_d*(sin(psi + psi_d*delta_t)-sin(psi)),
           v/psi_d*(-cos(psi + psi_d*delta_t)+cos(psi)),
           0,
           psi_d*delta_t,
           0;

      //std::cout << "F = " << std::endl << F << std::endl;
    }
    else
    {
      F << v*cos(psi)*delta_t,
           v*sin(psi)*delta_t,
           0,
           psi_d*delta_t,
           0;
      //std::cout << "F = " << std::endl << F << std::endl;
    }

    VectorXd n(n_x_);
    n.fill(0.0);
    n << 0.5*dt_2*cos(psi)*nu_a,
         0.5*dt_2*sin(psi)*nu_a,
         delta_t*nu_a,
         0.5*dt_2*nu_psidd,
         delta_t*nu_psidd;

    //cout << "n = " << endl << n << endl;


    Xsig_pred_.col(i) = x + F + n;

    while (Xsig_pred_.col(i)[3]> M_PI) Xsig_pred_.col(i)[3] -= 2.*M_PI;
    while (Xsig_pred_.col(i)[3]<-M_PI) Xsig_pred_.col(i)[3] += 2.*M_PI;


  }

  //cout << "Xsig_pred = " << endl << Xsig_pred_ <<endl;

  PredictMeanAndCovariance(&x_, &P_);
  //cout << x_ << endl;
  //cout << P_ << endl;
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
  int n_z = 2;
  VectorXd z = meas_package.raw_measurements_;
  //cout << "z = " << endl << z << endl;

  MatrixXd S(n_z, n_z);
  S.fill(0.0);

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  
  MatrixXd Zsig(n_z, 2*n_aug_ + 1);
  Zsig.fill(0.0);

  PredictLaserMeasurements(&z_pred, &S, &Zsig);

  NIS_laser_ = (z - z_pred).transpose()*S.inverse()*(z - z_pred);

  //if (NIS_laser_ > 5.991)
  //  return;

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //cout << z_diff << endl;
    //angle normalization
    //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  //cout << "Tc = " << endl << Tc << endl;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  //cout << "K = " << endl << K << endl;

  //residual
  VectorXd z_diff = z - z_pred;
  //cout << "z_diff = " << endl << z_diff << endl;

  //angle normalization
  //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //cout << "x_" << endl << x_ << endl;
  //cout << "P_" << endl << P_ << endl;

  //cout << "NIS_Laser_ = " << NIS_laser_ << endl;
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
  int n_z = 3;

  VectorXd z = meas_package.raw_measurements_;

  MatrixXd Zsig(n_z, 2*n_aug_ + 1);
  Zsig.fill(0.0);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);


  PredictRadarMeasurement(&z_pred, &S, &Zsig);

  this->NIS_radar_ = (z - z_pred).transpose()*S.inverse()*(z - z_pred);

  //if (NIS_radar_ > 7.815)
  //  return;

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //cout << z_diff << endl;
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
  //cout << "Tc = " << endl << Tc << endl;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  //cout << "K = " << endl << K << endl;

  //residual
  VectorXd z_diff = z - z_pred;
  //cout << "z_diff = " << endl << z_diff << endl;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  //cout << "x_" << endl << x_ << endl;
  //cout << "P_" << endl << P_ << endl;


  
  //cout << "NIS_radar_ = " << NIS_radar_ << endl;
}

void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out)
{
  //Augmenting the mean vector
  VectorXd x_aug(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  x_aug[5] = 0;//this->std_a_;
  x_aug[6] = 0;//this->std_yawdd_;

  //Augmenting the covariance matrix
  MatrixXd Q(n_aug_ - n_x_, n_aug_ - n_x_);
  Q << std_a_*std_a_,          0,
             0,       std_yawdd_*std_yawdd_;

  MatrixXd P_aug(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  //cout << P_aug << endl;
  P_aug.bottomRightCorner(n_aug_ - n_x_, n_aug_ - n_x_) = Q;
  
  //cout << "P_aug = " << endl << P_aug << endl;

  MatrixXd L = P_aug.llt().matrixL();
  
  //cout << "L = " << endl << L << endl;

  

  MatrixXd Xsig_aug(n_aug_, 2*n_aug_ + 1);

  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_)*L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*L.col(i);
  }

  //cout << "Xsig_aug = " << endl << Xsig_aug << endl;

  *Xsig_out = Xsig_aug;

}

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out)
{
  VectorXd x(n_x_);
  x.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++)
  {
    x += weights_[i]*Xsig_pred_.col(i);
  }

  while (x(3)> M_PI) x(3)-=2.*M_PI;
  while (x(3)<-M_PI) x(3)+=2.*M_PI;

  MatrixXd P(n_x_, n_x_);
  P.fill(0.0);
  for(int i = 0; i < 2*n_aug_ + 1; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //cout << "Xsig_pred_" << endl << Xsig_pred_ << endl;
    //cout << "x_diff(3) = " << x_diff(3) << endl;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P +=weights_[i]*(x_diff)*(x_diff).transpose();
  }

  //cout << "X_aug predicted = " << endl << x << endl;
  //cout << "P_aug predicted = " << endl << P << endl;

  *x_out = x;
  *P_out = P;
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out)
{
  //radar measurement dimensions
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    if(abs(Zsig(0,i)) > 0.001)
      Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    else
      Zsig(2,i) = 0;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;


  //print result
  //std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  //std::cout << "S: " << std::endl << S << std::endl;
  //std::cout << "Zsig: " << std::endl << Zsig << std::endl;

  //write result
  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}

void UKF::PredictLaserMeasurements(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out)
{
  int n_z = 2;

  VectorXd z_pred(n_z);
  z_pred.fill(0.0);

  MatrixXd S(n_z, n_z);
  S.fill(0.0);

  MatrixXd Zsig(n_z, 2*n_aug_ + 1);
  Zsig.fill(0.0);

  for(int i = 0; i < 2*n_aug_ + 1; i++)
  {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);

    Zsig(0,i) = px;
    Zsig(1,i) = py;
  }

  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_;
  S = S + R;


  //print result
  //std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  //std::cout << "S: " << std::endl << S << std::endl;
  //std::cout << "Zsig: " << std::endl << Zsig << std::endl;

  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}

void UKF::SetAccNoise(double std)
{
  this->std_a_ = std;
}

void UKF::SetPhiDDNoise(double std)
{
  this->std_yawdd_ = std;
}