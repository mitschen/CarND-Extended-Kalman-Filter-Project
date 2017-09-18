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
}
