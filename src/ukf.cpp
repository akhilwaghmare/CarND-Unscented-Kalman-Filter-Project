#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.setIdentity();

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.55;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  Xsig_pred_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  n_x_ = x_.rows();

  n_aug_ = n_x_ + 2;

  lambda_ = 3 - n_aug_;

  n_sigma_ = 2 * n_aug_ + 1;

  weights_ = VectorXd(n_sigma_);
  weights_[0] = lambda_/(lambda_+n_aug_);
  for (int i=1; i<n_sigma_; i++) {  //2n+1 weights
	  double weight = 0.5/(n_aug_+lambda_);
	  weights_[i] = weight;
  }

  time_us_ = 0;
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

  if (!is_initialized_) {
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		float rho = meas_package.raw_measurements_[0];
		float phi = meas_package.raw_measurements_[1];
		x_[0] = rho*cos(phi);
		x_[1] = rho*sin(phi);
	} else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		x_[0] = meas_package.raw_measurements_[0];
		x_[1] = meas_package.raw_measurements_[1];
	}

	time_us_ = meas_package.timestamp_;
	is_initialized_ = true;
	return;
  }

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
	  UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
	  UpdateLidar(meas_package);
  }
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

  MatrixXd Xsig_aug = getAugmentedSigmaPoints();
  Xsig_pred_ = getSigmaPointPrediction(Xsig_aug, delta_t);
  getPredictedMeanAndCovariance();
}

MatrixXd UKF::getAugmentedSigmaPoints() {
	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	x_aug.setZero();
	x_aug.head(n_x_) = x_;

	P_aug.setZero();
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(5, 5) = std_a_*std_a_;
	P_aug(6, 6) = std_yawdd_*std_yawdd_;

	MatrixXd A = P_aug.llt().matrixL();

	Xsig_aug.col(0) = x_aug;
	float c = sqrt(lambda_ + n_aug_);
	for (int i=0; i<n_aug_; i++) {
		Xsig_aug.col(i+1) = x_aug + c*A.col(i);
		Xsig_aug.col(n_aug_+i+1) = x_aug - c*A.col(i);
	}

	return Xsig_aug;
}

MatrixXd UKF::getSigmaPointPrediction(MatrixXd Xsig_aug, double delta_t) {
	MatrixXd Xsig_pred = MatrixXd(n_x_, n_sigma_);
	VectorXd x = x_;

	for (int i=0; i < n_sigma_; i++) {
		VectorXd x = Xsig_aug.col(i);
		VectorXd nu_k = VectorXd(n_x_);
		nu_k << delta_t*delta_t*cos(x[3])*x[5] / 2,
		delta_t*delta_t*sin(x[3])*x[5] / 2,
		delta_t*x[5],
		delta_t*delta_t*x[6] / 2,
		delta_t*x[6];
		VectorXd v = VectorXd(n_x_);
		if (x[4] == 0) {
			v << x[2]*cos(x[3])*delta_t,
			x[2]*sin(x[3])*delta_t,
			0,
			x[4]*delta_t,
			0;
		} else {
			v << (sin(x[3] + delta_t*x[4]) - sin(x[3])) * x[2] / x[4],
			(- cos(x[3] + delta_t*x[4]) + cos(x[3])) * x[2] / x[4],
			0,
			x[4]*delta_t,
			0;
		}

		VectorXd x_k1 = x.head(n_x_) + v + nu_k;
		Xsig_pred.col(i) = x_k1;
	}

	return Xsig_pred;
}

void UKF::getPredictedMeanAndCovariance() {
	x_.setZero();
	for (int i=0; i<n_sigma_; i++) {
		x_ = x_ + weights_[i]*Xsig_pred_.col(i);
	}

	P_.setZero();
	for (int i=0; i<n_sigma_; i++) {
		VectorXd y = Xsig_pred_.col(i) - x_;
		while (y(3)> M_PI) y(3)-=2.*M_PI;
		while (y(3)<-M_PI) y(3)+=2.*M_PI;
		P_ = P_ + weights_[i]*y*y.transpose();
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

  int n_z = 2;

  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0],
	meas_package.raw_measurements_[1];

  MatrixXd Zsig = predictSigmaPointsInLaserSpace();
  VectorXd z_pred = getPredictedMeanForLaser(Zsig);
  MatrixXd S = getPredictedCovarianceForLaser(Zsig, z_pred);

  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc.setZero();
  for (int i=0; i<n_sigma_; i++) {
	  VectorXd z_diff = Zsig.col(i) - z_pred;

	  VectorXd x_diff = Xsig_pred_.col(i) - x_;
	  while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
	  while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

	  Tc = Tc + weights_[i]*x_diff*z_diff.transpose();
  }

  MatrixXd K = Tc*S.inverse();

  x_ = x_ + K*(z - z_pred);
  P_ = P_ - K*S*K.transpose();
}

MatrixXd UKF::predictSigmaPointsInLaserSpace() {
	int n_z = 2;

	MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

	for(int i=0; i<n_sigma_; i++) {
		VectorXd x = Xsig_pred_.col(i);
		VectorXd measurement = VectorXd(n_z);
		measurement << x[0],
					   x[1];
		Zsig.col(i) = measurement;
	}

	return Zsig;
}

VectorXd UKF::getPredictedMeanForLaser(MatrixXd Zsig) {
	int n_z = 2;

	VectorXd z_pred = VectorXd(n_z);

	z_pred.setZero();
	for (int i=0; i<n_sigma_; i++) {
		z_pred = z_pred + weights_[i]*Zsig.col(i);
	}

	return z_pred;
}

MatrixXd UKF::getPredictedCovarianceForLaser(MatrixXd Zsig, VectorXd z_pred) {
	int n_z = 2;

	MatrixXd S = MatrixXd(n_z,n_z);

	MatrixXd R = MatrixXd(2,2);
	R.setZero();
	R(0,0) = std_laspx_*std_laspx_;
	R(1,1) = std_laspy_*std_laspy_;

	S.setZero();
	for (int i=0; i<n_sigma_; i++) {
		VectorXd y = Zsig.col(i) - z_pred;
		S = S + weights_[i]*y*y.transpose();
	}
	S = S + R;

	return S;
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

  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0],
	   meas_package.raw_measurements_[1],
	   meas_package.raw_measurements_[2];

  MatrixXd Zsig = predictSigmaPointsInRadarSpace();
  VectorXd z_pred = getPredictedMeanForRadar(Zsig);
  MatrixXd S = getPredictedCovarianceForRadar(Zsig, z_pred);

  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc.setZero();
  for (int i=0; i<n_sigma_; i++) {
	  VectorXd z_diff = Zsig.col(i) - z_pred;
	  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	  VectorXd x_diff = Xsig_pred_.col(i) - x_;
	  while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
	  while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

	  Tc = Tc + weights_[i]*x_diff*z_diff.transpose();
  }

  MatrixXd K = Tc*S.inverse();

  x_ = x_ + K*(z - z_pred);
  P_ = P_ - K*S*K.transpose();
}

MatrixXd UKF::predictSigmaPointsInRadarSpace() {
	int n_z = 3;

	MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

	for(int i=0; i<n_sigma_; i++) {
		VectorXd x = Xsig_pred_.col(i);
		VectorXd measurement = VectorXd(n_z);
		measurement << sqrt(x[0]*x[0] + x[1]*x[1]),
		atan2(x[1],x[0]),
		(x[0]*cos(x[3])*x[2] + x[1]*sin(x[3])*x[2]) / (sqrt(x[0]*x[0] + x[1]*x[1]));
		Zsig.col(i) = measurement;
	}

	return Zsig;
}

VectorXd UKF::getPredictedMeanForRadar(MatrixXd Zsig) {
	int n_z = 3;

	VectorXd z_pred = VectorXd(n_z);

	z_pred.setZero();
	for (int i=0; i<n_sigma_; i++) {
		z_pred = z_pred + weights_[i]*Zsig.col(i);
	}

	return z_pred;
}

MatrixXd UKF::getPredictedCovarianceForRadar(MatrixXd Zsig, VectorXd z_pred) {
	int n_z = 3;

	MatrixXd S = MatrixXd(n_z,n_z);

	MatrixXd R = MatrixXd(3,3);
	R.setZero();
	R(0,0) = std_radr_*std_radr_;
	R(1,1) = std_radphi_*std_radphi_;
	R(2,2) = std_radrd_*std_radrd_;

	S.setZero();
	for (int i=0; i<n_sigma_; i++) {
		VectorXd y = Zsig.col(i) - z_pred;
		while (y(1)> M_PI) y(1)-=2.*M_PI;
		while (y(1)<-M_PI) y(1)+=2.*M_PI;
		S = S + weights_[i]*y*y.transpose();
	}
	S = S + R;

	return S;
}
