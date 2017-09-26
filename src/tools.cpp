#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  int const c_maxCnt(estimations.size()>ground_truth.size()?
        estimations.size():ground_truth.size());
  assert(ground_truth.size()!=0);

  //set the number of RMSE elements expected
  VectorXd RMSE_return(ground_truth[0].size());
  RMSE_return.setZero();
  //calculate the sum of the squares
  for(int idx(0); idx < c_maxCnt; idx++)
  {
    VectorXd square_balance(estimations[idx]-ground_truth[idx]);
    square_balance = square_balance.array()*square_balance.array();
    RMSE_return = RMSE_return + square_balance;
  }
  //devide the sum with the number of elements
  RMSE_return = RMSE_return / c_maxCnt;
  //return the square-root
//  std::cout << "RMSE is "<<std::endl<<RMSE_return<<std::endl<<"**********"<<std::endl;

  return RMSE_return.array().sqrt();
}
