/*
 * PoseEKF.h
 *
 *  Created on: Oct 31, 2017
 *      Author: kevin
 *
 */
#ifndef INVIO_INCLUDE_INVIO_POSEEKF_H_
#define INVIO_INCLUDE_INVIO_POSEEKF_H_

#include <ros/ros.h>
#include <sophus/se3.hpp>
#include "../invio/vioParams.h"
#include "sensor_msgs/Imu.h"

// x = x, y, z, qw, qx, qy, qz, lambda, b_dx, b_dy, b_dz, b_wx, b_wy, b_wz, b_ax, b_ay, b_az, gx, gy, gz, bias_accel_x, bias_accel_y, bias_accel_z, bias_gyro_x, bias_gyro_y, bias_gyro_z;
#define STATE_SIZE 26

class PoseEKF {

	//this EKF allows invio to fuse in an imu to the motion estimate
	//NOTE: if running with dataset must be using sim time to initialize the ekf properly

public:
	PoseEKF();
	PoseEKF(ros::Time start);
	virtual ~PoseEKF();

	struct State{
		Eigen::Matrix<double, STATE_SIZE, 1> x; // the state
		ros::Time t; // current time of state estimate
		Eigen::Matrix<double, STATE_SIZE, STATE_SIZE> Sigma; // covariance

		/*
		 * returns unscaled position
		 */
		Eigen::Vector3d getPosition(){return Eigen::Vector3d(x(0), x(1), x(2));}
		Eigen::Quaterniond getQuat(){return Eigen::Quaterniond(x(3), x(4), x(5), x(6));}
		double getLambda(){return x(7);}
		Eigen::Vector3d getBodyFrameVelocity(){return Eigen::Vector3d(x(8), x(9), x(10));}
		Eigen::Vector3d getOmega(){return Eigen::Vector3d(x(11), x(12), x(13));}
		Eigen::Vector3d getBodyFrameAcceleration(){return Eigen::Vector3d(x(14), x(15), x(16));}
		Eigen::Vector3d getGravityVector(){return Eigen::Vector3d(x(17), x(18), x(19));}
		Eigen::Vector3d getAccelerometerBiases(){return Eigen::Vector3d(x(20), x(21), x(22));}
		Eigen::Vector3d getGyroscopeBiases(){return Eigen::Vector3d(x(23), x(24), x(25));}

		void setPosition(Eigen::Vector3d in){x(0)=in.x(); x(1)=in.y(); x(2)=in.z();}
		void setQuat(Eigen::Quaterniond in){x(3)=in.w(); x(4)=in.x(); x(5)=in.y(); x(6)=in.z();}
		void setLambda(double in){x(7)=in;}
		void setBodyFrameVelocity(Eigen::Vector3d in){x(8)=in.x(); x(9)=in.y(); x(10)=in.z();}
		void setOmega(Eigen::Vector3d in){x(11)=in.x(); x(12)=in.y(); x(13)=in.z();}
		void setBodyFrameAcceleration(Eigen::Vector3d in){x(14)=in.x(); x(15)=in.y(); x(16)=in.z();}
		void setGravityVector(Eigen::Vector3d in){x(17)=in.x(); x(18)=in.y(); x(19)=in.z();}
		void setAccelerometerBias(Eigen::Vector3d in){x(20)=in.x(); x(21)=in.y(); x(22)=in.z();}
		void setGyroscopeBias(Eigen::Vector3d in){x(23)=in.x(); x(24)=in.y(); x(25)=in.z();}


	} state;

	void process(ros::Time to_time); // predict the state to this time

	void updateWithVOPose(Sophus::SE3d pose, Eigen::Matrix<double, 6, 6> cov, ros::Time t_measured); // update with a pose estimate where cov is in order: [x, y, z, r, p, yaw]

	void updateWithIMU(sensor_msgs::Imu msg);

private:

	void resetState(); // place the state and covariance at initial conditions

	void computeStateTransitionJacobian(State& from_state, double dt, Eigen::Matrix<double, STATE_SIZE, STATE_SIZE>& F);

};

#endif /* INVIO_INCLUDE_INVIO_POSEEKF_H_ */
