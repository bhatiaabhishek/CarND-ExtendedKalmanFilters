#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include "math.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
	      0, 1, 0, 0;
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    x, P, F, H, R, Q
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
   ekf_.F_ = MatrixXd(4, 4);
   ekf_.F_ << 1, 0, 1, 0,
	      0, 1, 0, 1,
	      0, 0, 1, 0,
	      0, 0, 0, 1;
   ekf_.Q_ = MatrixXd(4,4);
   //state covariance matrix P
   ekf_.P_ = MatrixXd(4, 4);
   ekf_.P_ << 1, 0, 0, 0,
	      0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;
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
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float ro_v = measurement_pack.raw_measurements_(2);
      ekf_.x_ << ro*cos(phi), ro*sin(phi), 0,0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

     if (fabs(ekf_.x_(0)) < 0.001 and fabs(ekf_.x_(1)) < 0.001){
		ekf_.x_(0) = 0.001;
		ekf_.x_(1) = 0.001;
	}
    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
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
  //compute the time elapsed between the current and previous measurements

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.F_ << 1, 0, dt, 0,
	     0, 1, 0, dt,
	     0, 0, 1, 0,
	     0, 0, 0, 1;
	
  float dt_2 = pow(dt,2);
  float dt_3 = pow(dt,3);
  float dt_4 = pow(dt,4);
  float noise_ax_2 = 9;
  float noise_ay_2 = 9;
  ekf_.Q_ << dt_4*noise_ax_2/4, 0, dt_3*noise_ax_2/2, 0,
	          0, dt_4*noise_ay_2/4, 0, dt_3*noise_ay_2/2,
	          dt_3*noise_ax_2/2, 0, dt_2*noise_ax_2, 0,
	          0, dt_3*noise_ay_2/2, 0, dt_2*noise_ay_2;
  ekf_.Predict();
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  VectorXd M;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
        //recover state parameters
	float px = ekf_.x_(0);
	float py = ekf_.x_(1);
	float vx = ekf_.x_(2);
	float vy = ekf_.x_(3);


	//check division by zero
	float rho = sqrt(pow(px,2) + pow(py,2));
	float rho_2 = pow(rho,2);
	if (rho == 0) {
	    cout << "Division by zero";
	}
	else {
            if (fabs(px) < 0.001) { 
              if (px > 0) { px = 0.001;}
              else { px = -0.001; }
               }
            if (fabs(py) < 0.001) { 
              if (py > 0) { py = 0.001;}
              else { py = -0.001; }
               }
	    Hj_  << px/rho, py/rho, 0, 0,
	           -py/rho_2, px/rho_2, 0, 0,
	           py*((vx*py) - (vy*px))/pow(rho,3), px*((vy*px) - (vx*py))/pow(rho,3), px/rho, py/rho;
        ekf_.h_x = VectorXd(3);
        ekf_.h_x << rho, atan2(py,px), ((px*vx) + (py*vy))/rho;
        ekf_.H_ = Hj_;
        ekf_.R_ = R_radar_;
        M = VectorXd(3);
        M << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1],measurement_pack.raw_measurements_[2];
        ekf_.UpdateEKF(M);
	}
  } else {
    // Laser updates
    M = VectorXd(2);
    M << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
    ekf_.H_ = H_laser_; 
    ekf_.R_ = R_laser_;
    ekf_.Update(M);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
