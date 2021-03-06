#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
	rmse << 0,0,0,0;

  if(estimations.size() == 0 || ground_truth.size() == 0 
    || estimations.size() != ground_truth.size()) {
      cout << "Calculating RMSE, something went wrong";
      return rmse;
  }
  for (int i = 0; i < estimations.size(); i++) {
    VectorXd aux = estimations[i] - ground_truth[i];
    aux = aux.array()*aux.array();
    rmse += aux;
  }

  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
	// if (rmse[0] > 0.11 || rmse[1] > 0.11 || rmse[2] > 0.52 || rmse[3] > 0.52){
	// 	std::cout << "RMSE TOO HIGH" << endl;
	// }
  return rmse;
}