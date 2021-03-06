#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter()
  : x_(4)
  , P_(/*tbd*/)
  , F_(/*tbd*/)
  , Q_(/*tbd*/)
  , H_(/*tbd*/)
  , R_(/*tbd*/)
  {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  //estimation of the state, ignore the noise (mean = 0)
  x_ = F_ * x_;
  //estimate the covariance of the estimation (considering uncertainty of the noise)
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  MatrixXd const Ht = H_.transpose();

  //calculate the error of prediction and measurement
  VectorXd const y = z - H_ * x_;
  MatrixXd const S = H_ * P_ * Ht + R_;

  //calculate the kalman gain
  MatrixXd const K = P_ * Ht * S.inverse();

  //calculate the final state and covariance
  x_ = x_ + K * y;
  P_ = (MatrixXd::Identity(K.rows(),H_.cols()) - K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {


  //translate predicted state x' into the measurement space
  double const pred_rho(sqrt(x_(0)*x_(0) + x_(1)*x_(1)));
  double const pred_phi(atan2(x_(1),x_(0)));

  //sanity check
  bool const avoidZeroDivision(fabs(pred_rho)<0.0001);
  double const pred_rhodot(avoidZeroDivision?0.: (x_(0)*x_(2)+x_(1)*x_(3)) / pred_rho);


  //summarize the prediction vector
  VectorXd prediction(3);
  prediction << pred_rho, pred_phi, pred_rhodot;

  VectorXd y(z-prediction);
  const double pi(3.14159265);

  //adjust the range of phi
  y(1) = fmod(y(1),pi);

  //do filter update
  MatrixXd const Ht = H_.transpose();
  MatrixXd const S = H_ * P_ * Ht + R_;

  //calculate the kalman gain
  MatrixXd const K = P_ * Ht * S.inverse();

  //calculate the final state and covariance
  x_ = x_ + K * y;
  P_ = (MatrixXd::Identity(K.rows(),H_.cols()) - K*H_)*P_;

}
