#include <iostream>
#include "tools.h"
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

        VectorXd rmse(4);
	rmse << 0,0,0,0;

	VectorXd res(4);
	res << 0, 0, 0, 0;
	VectorXd err(4);
	float n;
	n = 0;
    if ((estimations.size() != 0) and (estimations.size() == ground_truth.size())) {
        
    
	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        // ... your code here
        n = n + 1;
        err = estimations[i] - ground_truth[i];
        err =  err.array()*err.array();
	    res = res + err;
	}
    
    
	//calculate the mean
	// ... your code here
    res = res.array()/n;
	//calculate the squared root
	// ... your code here
    rmse = res.array().sqrt();
	//return the result
    }
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
        MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);


	//check division by zero
	float rho = sqrt(pow(px,2) + pow(py,2));
	float rho_2 = pow(rho,2);
	if (rho == 0) {
	    std::cout << "Division by zero";
	}
	else {
	    Hj  << px/rho, py/rho, 0, 0,
	           -py/rho_2, px/rho_2, 0, 0,
	           py*(vx*py - vy*px)/pow(rho,3), px*(vy*px - vx*py)/pow(rho,3), px/rho, py/rho;
	}
	//compute the Jacobian matrix

	return Hj;
}
