#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd RMSE_return(4);
  int const c_maxCnt(estimations.size()>ground_truth.size()?
        estimations.size():ground_truth.size());

  //calculate the sum of the squares
  for(int idx(0); idx < c_maxCnt; idx++)
  {
    VectorXd square_balance(estimations[idx]-ground_truth[idx]);
    square_balance = square_balance.array()*square_balance.array();
    RMSE_return = RMSE_return + square_balance;
  }
  //devide the sum with the number of elements
  RMSE_return = RMSE_return / c_maxCnt;
  //return the sqare-root
  return RMSE_return.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd H(3,4);
  double const x(x_state(0)), y(x_state(1))
            , vx(x_state(2)), vy(x_state(3));

  double const xy_2(x*x + y*y);
  if(fabs(xy_2) < 0.0001)
  {
    cout << "CalculateJacobian () - Error - Division by zero" << endl;
    return H;
  }
  double const v1(vx*y-vy*x);
  H(0,0) = x / sqrt(xy_2);
  H(0,1) = y / sqrt(xy_2);
  H(0,2) = 0.;
  H(0,3) = 0.;
  H(1,0) = -y / xy_2;
  H(1,1) = x / xy_2;
  H(1,2) = 0.;
  H(1,3) = 0.;
  H(2,0) = (y * v1)/ (sqrt(xy_2)*xy_2);
  H(2,1) = (-x * v1)/ (sqrt(xy_2)*xy_2);
  H(2,2) = x/sqrt(xy_2);
  H(2,3) = y/sqrt(xy_2);
  return H;
}
