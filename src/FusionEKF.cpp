#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <unistd.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF()
  : ekf_()
  , is_initialized_(false)
  , previous_timestamp_(0)
  , tools()
  , R_laser_(2, 2)
  , R_radar_(3, 3)
  , H_laser_(2, 4)
  , Hj_(3, 4)
{
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0     , 0.0225;
  //measurement matrix for laser
  //we consider only position measurement - no velocity given
  H_laser_ << 1, 0, 0, 0,
		      0, 1, 0, 0;
  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0     , 0,
              0   , 0.0009, 0,
              0   , 0     , 0.09;
  //measurement matrix for radar
  //Please note: this matrix must be calculated each iteration using Jacobian
  Hj_ << 0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/

  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    previous_timestamp_ = measurement_pack.timestamp_;
    //initial covariance state matrix
    //following suggestion of Chapter 11
    MatrixXd P(4,4);
    P << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1000, 0,
         0, 0, 0, 1000;

    //State transition matrix
    //For the first step, x' is simply x
    MatrixXd F(4, 4);
    F << 1, 0, 1, 0,
         0, 1, 0, 1,
         0, 0, 1, 0,
         0, 0, 0, 1;

    //Uncertainty Q
    MatrixXd Q(4,4);

    VectorXd X(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double const rho(measurement_pack.raw_measurements_[0]);
      double const phi(measurement_pack.raw_measurements_[1]);
      double const rhodot(measurement_pack.raw_measurements_[2]);
      double const x(cos(phi) * rho);
      double const y(sin(phi) * rho);
      X << x, y, 0, 0;
      ekf_.Init(X, P, F, Hj_, R_radar_, Q);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      double const x(measurement_pack.raw_measurements_[0]);
      double const y(measurement_pack.raw_measurements_[1]);
      //the initial state x vector contains only position, no velocity
      X << x, y, 0, 0;
      //initialize the KMF
      ekf_.Init(X, P, F, H_laser_, R_laser_, Q);
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  //identify the time difference
  double const dt( (measurement_pack.timestamp_ - previous_timestamp_)/ 1000000.0);
  double const dt2(dt*dt);
  double const dt3(dt2*dt);
  double const dt4(dt3*dt2);
  //update previous timestamp for next iteration
  previous_timestamp_ = measurement_pack.timestamp_;

  double const noise_ax(9.);
  double const noise_ay(9.);
  //calculate the process noise covariance matrix Q (each timestep, noise
  //		increases if timespan is longer)
  MatrixXd Q(4,4);
  Q << dt4/4*noise_ax, 0             , dt3/2*noise_ax, 0,
       0             , dt4/4*noise_ay, 0             , dt3/2*noise_ay,
       dt3/2*noise_ax, 0             , dt2*noise_ax  , 0,
       0             , dt3/2*noise_ay, 0             , dt2*noise_ay;
  //update the state transition matrix
  MatrixXd F(4,4);
  F << 1, 0, dt, 0,
       0, 1, 0, dt,
       0, 0, 1, 0,
       0, 0, 0, 1;
  //do prediction
  ekf_.F_ = F;
  ekf_.Q_ = Q;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    //update the state transition matrix and the covariance
    //use the prediction ekf_.x_ to calculate jacobian
    //Jacobian is making non-linearities linear again
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    //update the state transition matrix and the covariance
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
