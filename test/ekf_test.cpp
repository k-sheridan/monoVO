/*
 * ekf_test.cpp
 *
 *  Created on: Jan 1, 2018
 *      Author: kevin
 */


#include <ros/ros.h>

#include <PoseEKF.h>

int main(int argc, char **argv)
{
	ros::init(argc, argv, "ekf_test"); // initializes ros

	PoseEKF ekf;

	ekf.resetState();

	Eigen::Matrix<double, STATE_SIZE, STATE_SIZE> F;

	//-=-=-=
	ekf.computeStateTransitionJacobian(ekf.state, 0.1, F);
	ROS_INFO_STREAM("transition for starting state at dt = 0.1: " << F);

	//-=-=-=
	F.setZero();
	ekf.computeStateTransitionJacobian(ekf.state, 0.0, F);
	ROS_INFO_STREAM("transition for starting state at dt = 0.0: " << F);

	//-=-=-=
	ekf.state.setOmega(Eigen::Vector3d(1, 1, 1));
	F.setZero();
	ekf.computeStateTransitionJacobian(ekf.state, 0.1, F);
	ROS_INFO_STREAM("transition for non zero omega state at dt = 0.0: " << F);

	//-=-=-=-=-
	Eigen::Matrix<double, 7, 6> J;
	Eigen::Matrix<double, 6, 1> vec;
	vec << 0, 0, 0, 0, 0, 0;
	ekf.computeAngleAxis2QuaternionSpaceJacobian(vec, J);
	ROS_INFO_STREAM("test small angle a2q jacobain: " << J);

	//-=-=-=-=-
	J.setZero();
	vec << 0, 0, 0, 0.1, 0.2, 0.3;
	ekf.computeAngleAxis2QuaternionSpaceJacobian(vec, J);
	ROS_INFO_STREAM("test normal angle a2q jacobain: " << J);

	return 0;
}

