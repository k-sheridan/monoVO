/*
 * PoseEKF.cpp
 *
 *  Created on: Oct 31, 2017
 *      Author: kevin
 */

#include <PoseEKF.h>

PoseEKF::PoseEKF() {
	// TODO Auto-generated constructor stub
	this->resetState();
}

PoseEKF::PoseEKF(ros::Time start)
{
	this->state.t = start;
	this->resetState();
}

PoseEKF::~PoseEKF() {
	// TODO Auto-generated destructor stub
}

/*
 * set the initial state and initial uncertainties
 */
void PoseEKF::resetState()
{

}

/*
 * this predicts the state forward in time and updates the uncertainties and correlations due to prediction uncertainty
 */
void PoseEKF::process(ros::Time to_time)
{
	double dt = (to_time - this->state.t).toSec();
	ROS_ASSERT(dt >= 0);

	// first compute the linearized transition function at the previous state
	Eigen::Matrix<double, STATE_SIZE, STATE_SIZE> F;
	//TODO check if F will be initialized with all zeros
	this->computeStateTransitionJacobian(this->state, dt, F);

	// now we can predict the state forward with the nonlinear update function

	//r = r + R(q)*((lambda)*(dt*b_r' + 0.5*dt^2*b_r'')) { this is where the scale is estimated }
	this->state.setPosition(this->state.getPosition() + this->state.getQuat() * (this->state.getLambda() * (dt * this->state.getBodyFrameVelocity() + 0.5*dt*dt * this->state.getBodyFrameAcceleration())));

	//q = q*d_q(b_w*dt)
	// and
	//g = d_q(b_w*dt)*g
	double omega_norm = this->state.getOmega().norm();
	double theta = dt * omega_norm;
	if(omega_norm > 1e-10){ // this ensures that we don't divide by zero. im omega is zero that means that the delta quaternion is the identity quaternion
		Eigen::Vector3d v = this->state.getOmega() / omega_norm;

		double sine_theta_2 = sin(theta/2.0);

		Eigen::Quaterniond dq = Eigen::Quaterniond(cos(theta/2.0), v.x()*sine_theta_2, v.y()*sine_theta_2, v.z()*sine_theta_2);

		this->state.setQuat(this->state.getQuat() * dq); // finally rotate the quat by the delta quat
		this->state.setGravityVector(dq.inverse() * this->state.getGravityVector()); // multiply the gravity vector by the conjugate of the quaternion

		ROS_ASSERT(this->state.getQuat().norm() == 1.0);
	}

	//lambda^(-1) = lambda^(-1)

	//b_r' = b_r' + dt*b_dr''
	this->state.setBodyFrameVelocity(this->state.getBodyFrameVelocity() + (dt * this->state.getBodyFrameAcceleration()));

	//b_w = b_w

	//b_r'' = b_r''

	//biases = biases
}

/*
 * need to convert the euler covariance matrix to a quaternion covariance matrix before update can
 * be performed. This is a completely linear process other wise
 */
void PoseEKF::updateWithVOPose(Sophus::SE3d pose, Eigen::Matrix<double, 6, 6> cov, ros::Time t_measured)
{

}

/*
 * due to the scale factor, and quaternion, and body frame velocity... this is not a linear process
 */
void PoseEKF::computeStateTransitionJacobian(State& from_state, double dt, Eigen::Matrix<double, STATE_SIZE, STATE_SIZE>& F)
{

}

