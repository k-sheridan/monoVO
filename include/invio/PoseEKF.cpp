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

PoseEKF::PoseEKF(ros::Time start) {
	this->state.t = start;
	this->resetState();
}

PoseEKF::~PoseEKF() {
	// TODO Auto-generated destructor stub
}

/*
 * set the initial state and initial uncertainties
 */
void PoseEKF::resetState() {

}

/*
 * this predicts the state forward in time and updates the uncertainties and correlations due to prediction uncertainty
 */
void PoseEKF::process(ros::Time to_time) {
	double dt = (to_time - this->state.t).toSec();
	ROS_ASSERT(dt >= 0);

	// first compute the linearized transition function at the previous state
	Eigen::Matrix<double, STATE_SIZE, STATE_SIZE> F;
	F.setZero(); // ensure that this matrix is zeros
	//TODO check if F will be initialized with all zeros
	this->computeStateTransitionJacobian(this->state, dt, F);

	// now we can predict the state forward with the nonlinear update function

	//r = r + R(q)*((lambda)*(dt*b_r' + 0.5*dt^2*b_r'')) { this is where the scale is estimated }
	this->state.setPosition(
			this->state.getPosition()
					+ this->state.getQuat()
							* (this->state.getLambda()
									* (dt * this->state.getBodyFrameVelocity()
											+ 0.5 * dt * dt
													* this->state.getBodyFrameAcceleration())));

	//q = q*d_q(b_w*dt)
	// and
	//g = d_q(b_w*dt)*g
	double omega_norm = this->state.getOmega().norm();
	double theta = dt * omega_norm;
	if (omega_norm > 1e-10) { // this ensures that we don't divide by zero. im omega is zero that means that the delta quaternion is the identity quaternion
		Eigen::Vector3d v = this->state.getOmega() / omega_norm;

		double sine_theta_2 = sin(theta / 2.0);

		Eigen::Quaterniond dq = Eigen::Quaterniond(cos(theta / 2.0),
				v.x() * sine_theta_2, v.y() * sine_theta_2,
				v.z() * sine_theta_2);

		this->state.setQuat(this->state.getQuat() * dq); // finally rotate the quat by the delta quat
		this->state.setGravityVector(
				dq.inverse() * this->state.getGravityVector()); // multiply the gravity vector by the conjugate of the quaternion

		ROS_ASSERT(this->state.getQuat().norm() == 1.0);
	}

	//lambda^(-1) = lambda^(-1)

	//b_r' = b_r' + dt*b_dr''
	this->state.setBodyFrameVelocity(
			this->state.getBodyFrameVelocity()
					+ (dt * this->state.getBodyFrameAcceleration()));

	//b_w = b_w

	//b_r'' = b_r''

	//biases = biases
}

/*
 * need to convert the euler covariance matrix to a quaternion covariance matrix before update can
 * be performed. This is a completely linear process other wise
 */
void PoseEKF::updateWithVOPose(Sophus::SE3d pose,
		Eigen::Matrix<double, 6, 6> cov, ros::Time t_measured) {

}

/*
 * due to the scale factor, and quaternion, and body frame velocity... this is not a linear process
 */
void PoseEKF::computeStateTransitionJacobian(State& from_state, double dt,
		Eigen::Matrix<double, STATE_SIZE, STATE_SIZE>& F) {

	double x = from_state.x(0);
	double y = from_state.x(1);
	double z = from_state.x(2);
	double qw = from_state.x(3);
	double qx = from_state.x(4);
	double qy = from_state.x(5);
	double qz = from_state.x(6);
	double lambda = from_state.x(7);
	double b_dx = from_state.x(8);
	double b_dy = from_state.x(9);
	double b_dz = from_state.x(10);
	double b_wx = from_state.x(11);
	double b_wy = from_state.x(12);
	double b_wz = from_state.x(13);
	double b_ax = from_state.x(14);
	double b_ay = from_state.x(15);
	double b_az = from_state.x(16);
	double gx = from_state.x(17);
	double gy = from_state.x(18);
	double gz = from_state.x(19);

	double t2 = b_dy * dt;
	double t3 = dt * dt;
	double t4 = b_ay * t3 * (1.0 / 2.0);
	double t5 = t2 + t4;
	double t6 = b_az * t3 * (1.0 / 2.0);
	double t7 = t2 + t6;
	double t8 = b_dx * dt;
	double t9 = b_ax * t3 * (1.0 / 2.0);
	double t10 = t8 + t9;
	double t11 = qy * qy;
	double t12 = t11 * 2.0;
	double t13 = qz * qz;
	double t14 = t13 * 2.0;
	double t15 = t12 + t14 - 1.0;
	double t16 = qw * qy * 2.0;
	double t21 = qx * qz * 2.0;
	double t17 = t16 - t21;
	double t18 = qw * qz * 2.0;
	double t19 = qx * qy * 2.0;
	double t20 = t18 + t19;
	double t22 = lambda * qx * t7 * 2.0;
	double t23 = lambda * qz * t7 * 2.0;
	double t24 = t18 - t19;
	double t25 = qw * qx * 2.0;
	double t26 = qy * qz * 2.0;
	double t27 = t25 + t26;
	double t28 = qx * qx;
	double t29 = t28 * 2.0;
	double t30 = t14 + t29 - 1.0;
	double t31 = lambda * qx * t5 * 2.0;
	double t32 = lambda * qy * t10 * 2.0;
	double t33 = lambda * qw * t5 * 2.0;
	double t34 = lambda * qz * t5 * 2.0;
	double t35 = lambda * qx * t10 * 2.0;
	double t36 = lambda * qy * t5 * 2.0;
	double t37 = t16 + t21;
	double t38 = t25 - t26;
	double t39 = t12 + t29 - 1.0;
	double t40 = fabs(b_wx);
	double t41 = fabs(b_wy);
	double t42 = fabs(b_wz);
	double t43 = t40 * t40;
	double t44 = t41 * t41;
	double t45 = t42 * t42;
	double t46 = t43 + t44 + t45;

	if (t46 != 0) {

		double t47 = sqrt(t46);
		double t48 = dt * t47 * (1.0 / 2.0);
		double t49 = sin(t48);
		double t50 = 1.0 / sqrt(t46);
		double t51 = (b_wx / fabs(b_wx));
		double t52 = 1.0 / pow(t46, 3.0 / 2.0);
		double t53 = cos(t48);
		double t54 = 1.0 / t46;
		double t55 = (b_wy / fabs(b_wy));
		double t56 = (b_wz / fabs(b_wz));
		double t57 = b_wz * t49 * t50;
		double t58 = b_wx * t49 * t50;
		double t59 = qw * t49 * t50;
		double t60 = b_wy * t49 * t50;
		double t61 = qy * t49 * t50;
		double t62 = t49 * t49;
		double t63 = 1.0 / (t46 * t46);
		double t64 = b_wy * b_wy;
		double t65 = b_wz * b_wz;
		double t66 = t53 * t53;
		double t67 = t49 * t50 * t53 * 2.0;
		double t68 = b_wx * t54 * t62 * 2.0;
		double t69 = t40 * t51 * t62 * t63 * t65 * 4.0;
		double t70 = b_wx * b_wx;
		double t71 = b_wy * t54 * t62 * 2.0;
		double t72 = b_wz * t40 * t49 * t51 * t52 * t53 * 2.0;
		double t73 = b_wz * dt * t40 * t51 * t54 * t62;
		double t74 = b_wx * b_wy * dt * t40 * t49 * t51 * t52 * t53 * 2.0;
		double t75 = dt * t41 * t49 * t52 * t53 * t55 * t65 * 2.0;
		double t76 = b_wz * t41 * t49 * t52 * t53 * t55 * 2.0;
		double t77 = b_wz * dt * t41 * t54 * t55 * t62;
		double t78 = b_wx * b_wy * dt * t41 * t49 * t52 * t53 * t55 * 2.0;
		double t79 = b_wz * t54 * t62 * 2.0;
		double t80 = b_wz * t54 * t62 * 4.0;
		double t81 = dt * t42 * t49 * t52 * t53 * t56 * t65 * 2.0;
		double t82 = b_wz * dt * t42 * t54 * t56 * t66;
		double t83 = b_wx * b_wy * t42 * t56 * t62 * t63 * 4.0;
		double t84 = b_wx * b_wy * t54 * t62 * 2.0;
		double t85 = b_wx * t54 * t62 * 4.0;
		double t86 = t40 * t51 * t62 * t63 * t64 * 4.0;
		double t87 = dt * t40 * t49 * t51 * t52 * t53 * t70 * 2.0;
		double t88 = b_wx * dt * t40 * t51 * t54 * t66;
		double t89 = b_wy * b_wz * t40 * t51 * t62 * t63 * 4.0;
		double t90 = b_wy * dt * t40 * t51 * t54 * t66;
		double t91 = b_wx * b_wz * dt * t40 * t49 * t51 * t52 * t53 * 2.0;
		double t92 = b_wy * t54 * t62 * 4.0;
		double t93 = dt * t41 * t49 * t52 * t53 * t55 * t70 * 2.0;
		double t94 = dt * t41 * t49 * t52 * t53 * t55 * t64 * 2.0;
		double t95 = b_wy * dt * t41 * t54 * t55 * t66;
		double t96 = b_wx * b_wz * dt * t41 * t49 * t52 * t53 * t55 * 2.0;
		double t97 = b_wx * t41 * t49 * t52 * t53 * t55 * 2.0;
		double t98 = b_wx * dt * t41 * t54 * t55 * t62;
		double t99 = b_wy * b_wz * dt * t41 * t49 * t52 * t53 * t55 * 2.0;
		double t100 = dt * t42 * t49 * t52 * t53 * t56 * t70 * 2.0;
		double t101 = dt * t42 * t49 * t52 * t53 * t56 * t64 * 2.0;
		double t102 = b_wy * dt * t42 * t54 * t56 * t66;
		double t103 = b_wx * b_wz * dt * t42 * t49 * t52 * t53 * t56 * 2.0;
		double t104 = b_wx * t42 * t49 * t52 * t53 * t56 * 2.0;
		double t105 = b_wx * dt * t42 * t54 * t56 * t62;
		double t106 = b_wy * b_wz * dt * t42 * t49 * t52 * t53 * t56 * 2.0;
		double t107 = b_wx * b_wz * t54 * t62 * 2.0;
		double t108 = b_wy * t49 * t50 * t53 * 2.0;
		double t109 = b_wy * b_wz * t54 * t62 * 2.0;
		F(0,0) = 1.0;
		F(0,3) = t34 - lambda * qy * t7 * 2.0;
		F(0,4) = t23 + t36;
		F(0,5) = t31 - lambda * qw * t7 * 2.0 - lambda * qy * t10 * 4.0;
		F(0,6) = t22 + t33 - lambda * qz * t10 * 4.0;
		F(0,7) = -t7 * t17 + t5 * t20 - t10 * t15;
		F(0,8) = -dt * lambda * t15;
		F(0,9) = -dt * lambda * t17 + dt * lambda * t20;
		F(0,14) = lambda * t3 * t15 * (-1.0 / 2.0);
		F(0,15) = lambda * t3 * t20 * (1.0 / 2.0);
		F(0,16) = lambda * t3 * t17 * (-1.0 / 2.0);
		F(1,1) = 1.0;
		F(1,3) = t22 - lambda * qz * t10 * 2.0;
		F(1,4) = t32 + lambda * qw * t7 * 2.0 - lambda * qx * t5 * 4.0;
		F(1,5) = t23 + t35;
		F(1,6) = lambda * qw * t10 * -2.0 + lambda * qy * t7 * 2.0
				- lambda * qz * t5 * 4.0;
		F(1,7) = t7 * t27 - t10 * t24 - t5 * t30;
		F(1,8) = -dt * lambda * t24;
		F(1,9) = dt * lambda * t27 - dt * lambda * t30;
		F(1,14) = lambda * t3 * t24 * (-1.0 / 2.0);
		F(1,15) = lambda * t3 * t30 * (-1.0 / 2.0);
		F(1,16) = lambda * t3 * t27 * (1.0 / 2.0);
		F(2,2) = 1.0;
		F(2,3) = -t31 + t32;
		F(2,4) = -t33 - lambda * qx * t7 * 4.0 + lambda * qz * t10 * 2.0;
		F(2,5) = t34 + lambda * qw * t10 * 2.0 - lambda * qy * t7 * 4.0;
		F(2,6) = t35 + t36;
		F(2,7) = -t5 * t38 - t7 * t39 + t10 * t37;
		F(2,8) = dt * lambda * t37;
		F(2,9) = -dt * lambda * t38 - dt * lambda * t39;
		F(2,14) = lambda * t3 * t37 * (1.0 / 2.0);
		F(2,15) = lambda * t3 * t38 * (-1.0 / 2.0);
		F(2,16) = lambda * t3 * t39 * (-1.0 / 2.0);
		F(3,3) = t53;
		F(3,4) = -b_wx * t49 * t50;
		F(3,5) = -b_wy * t49 * t50;
		F(3,6) = -b_wz * t49 * t50;
		F(3,11) = -qx * t49 * t50 + b_wx * qx * t40 * t49 * t51 * t52
				+ b_wy * qy * t40 * t49 * t51 * t52
				+ b_wz * qz * t40 * t49 * t51 * t52
				- dt * qw * t40 * t49 * t50 * t51 * (1.0 / 2.0)
				- b_wx * dt * qx * t40 * t51 * t53 * t54 * (1.0 / 2.0)
				- b_wy * dt * qy * t40 * t51 * t53 * t54 * (1.0 / 2.0)
				- b_wz * dt * qz * t40 * t51 * t53 * t54 * (1.0 / 2.0);
		F(3,12) = -qy * t49 * t50 + b_wx * qx * t41 * t49 * t52 * t55
				+ b_wy * qy * t41 * t49 * t52 * t55
				+ b_wz * qz * t41 * t49 * t52 * t55
				- dt * qw * t41 * t49 * t50 * t55 * (1.0 / 2.0)
				- b_wx * dt * qx * t41 * t53 * t54 * t55 * (1.0 / 2.0)
				- b_wy * dt * qy * t41 * t53 * t54 * t55 * (1.0 / 2.0)
				- b_wz * dt * qz * t41 * t53 * t54 * t55 * (1.0 / 2.0);
		F(3,13) = -qz * t49 * t50 + b_wx * qx * t42 * t49 * t52 * t56
				+ b_wy * qy * t42 * t49 * t52 * t56
				+ b_wz * qz * t42 * t49 * t52 * t56
				- dt * qw * t42 * t49 * t50 * t56 * (1.0 / 2.0)
				- b_wx * dt * qx * t42 * t53 * t54 * t56 * (1.0 / 2.0)
				- b_wy * dt * qy * t42 * t53 * t54 * t56 * (1.0 / 2.0)
				- b_wz * dt * qz * t42 * t53 * t54 * t56 * (1.0 / 2.0);
		F(4,3) = t58;
		F(4,4) = t53;
		F(4,5) = t57;
		F(4,6) = -b_wy * t49 * t50;
		F(4,11) = t59 - b_wx * qw * t40 * t49 * t51 * t52
				- b_wz * qy * t40 * t49 * t51 * t52
				+ b_wy * qz * t40 * t49 * t51 * t52
				- dt * qx * t40 * t49 * t50 * t51 * (1.0 / 2.0)
				+ b_wx * dt * qw * t40 * t51 * t53 * t54 * (1.0 / 2.0)
				+ b_wz * dt * qy * t40 * t51 * t53 * t54 * (1.0 / 2.0)
				- b_wy * dt * qz * t40 * t51 * t53 * t54 * (1.0 / 2.0);
		F(4,12) = -qz * t49 * t50 - b_wx * qw * t41 * t49 * t52 * t55
				- b_wz * qy * t41 * t49 * t52 * t55
				+ b_wy * qz * t41 * t49 * t52 * t55
				- dt * qx * t41 * t49 * t50 * t55 * (1.0 / 2.0)
				+ b_wx * dt * qw * t41 * t53 * t54 * t55 * (1.0 / 2.0)
				+ b_wz * dt * qy * t41 * t53 * t54 * t55 * (1.0 / 2.0)
				- b_wy * dt * qz * t41 * t53 * t54 * t55 * (1.0 / 2.0);
		F(4,13) = t61 - b_wx * qw * t42 * t49 * t52 * t56
				- b_wz * qy * t42 * t49 * t52 * t56
				+ b_wy * qz * t42 * t49 * t52 * t56
				- dt * qx * t42 * t49 * t50 * t56 * (1.0 / 2.0)
				+ b_wx * dt * qw * t42 * t53 * t54 * t56 * (1.0 / 2.0)
				+ b_wz * dt * qy * t42 * t53 * t54 * t56 * (1.0 / 2.0)
				- b_wy * dt * qz * t42 * t53 * t54 * t56 * (1.0 / 2.0);
		F(5,3) = t60;
		F(5,4) = -t57;
		F(5,5) = t53;
		F(5,6) = t58;
		F(5,11) = qz * t49 * t50 - b_wy * qw * t40 * t49 * t51 * t52
				+ b_wz * qx * t40 * t49 * t51 * t52
				- b_wx * qz * t40 * t49 * t51 * t52
				- dt * qy * t40 * t49 * t50 * t51 * (1.0 / 2.0)
				+ b_wy * dt * qw * t40 * t51 * t53 * t54 * (1.0 / 2.0)
				- b_wz * dt * qx * t40 * t51 * t53 * t54 * (1.0 / 2.0)
				+ b_wx * dt * qz * t40 * t51 * t53 * t54 * (1.0 / 2.0);
		F(5,12) = t59 - b_wy * qw * t41 * t49 * t52 * t55
				+ b_wz * qx * t41 * t49 * t52 * t55
				- b_wx * qz * t41 * t49 * t52 * t55
				- dt * qy * t41 * t49 * t50 * t55 * (1.0 / 2.0)
				+ b_wy * dt * qw * t41 * t53 * t54 * t55 * (1.0 / 2.0)
				- b_wz * dt * qx * t41 * t53 * t54 * t55 * (1.0 / 2.0)
				+ b_wx * dt * qz * t41 * t53 * t54 * t55 * (1.0 / 2.0);
		F(5,13) = -qx * t49 * t50 - b_wy * qw * t42 * t49 * t52 * t56
				+ b_wz * qx * t42 * t49 * t52 * t56
				- b_wx * qz * t42 * t49 * t52 * t56
				- dt * qy * t42 * t49 * t50 * t56 * (1.0 / 2.0)
				+ b_wy * dt * qw * t42 * t53 * t54 * t56 * (1.0 / 2.0)
				- b_wz * dt * qx * t42 * t53 * t54 * t56 * (1.0 / 2.0)
				+ b_wx * dt * qz * t42 * t53 * t54 * t56 * (1.0 / 2.0);
		F(6,3) = t57;
		F(6,4) = t60;
		F(6,5) = -t58;
		F(6,6) = t53;
		F(6,11) = -t61 - b_wz * qw * t40 * t49 * t51 * t52
				- b_wy * qx * t40 * t49 * t51 * t52
				+ b_wx * qy * t40 * t49 * t51 * t52
				- dt * qz * t40 * t49 * t50 * t51 * (1.0 / 2.0)
				+ b_wz * dt * qw * t40 * t51 * t53 * t54 * (1.0 / 2.0)
				+ b_wy * dt * qx * t40 * t51 * t53 * t54 * (1.0 / 2.0)
				- b_wx * dt * qy * t40 * t51 * t53 * t54 * (1.0 / 2.0);
		F(6,12) = qx * t49 * t50 - b_wz * qw * t41 * t49 * t52 * t55
				- b_wy * qx * t41 * t49 * t52 * t55
				+ b_wx * qy * t41 * t49 * t52 * t55
				- dt * qz * t41 * t49 * t50 * t55 * (1.0 / 2.0)
				+ b_wz * dt * qw * t41 * t53 * t54 * t55 * (1.0 / 2.0)
				+ b_wy * dt * qx * t41 * t53 * t54 * t55 * (1.0 / 2.0)
				- b_wx * dt * qy * t41 * t53 * t54 * t55 * (1.0 / 2.0);
		F(6,13) = t59 - b_wz * qw * t42 * t49 * t52 * t56
				- b_wy * qx * t42 * t49 * t52 * t56
				+ b_wx * qy * t42 * t49 * t52 * t56
				- dt * qz * t42 * t49 * t50 * t56 * (1.0 / 2.0)
				+ b_wz * dt * qw * t42 * t53 * t54 * t56 * (1.0 / 2.0)
				+ b_wy * dt * qx * t42 * t53 * t54 * t56 * (1.0 / 2.0)
				- b_wx * dt * qy * t42 * t53 * t54 * t56 * (1.0 / 2.0);
		F(7,7) = 1.0;
		F(8,8) = 1.0;
		F(8,14) = dt;
		F(9,9) = 1.0;
		F(9,15) = dt;
		F(10,9) = 1.0;
		F(10,16) = dt;
		F(11,11) = 1.0;
		F(12,12) = 1.0;
		F(13,13) = 1.0;
		F(14,14) = 1.0;
		F(15,15) = 1.0;
		F(16,16) = 1.0;
		F(17,11) = gx
				* (t69 + t86 - dt * t40 * t49 * t51 * t52 * t53 * t64 * 2.0
						- dt * t40 * t49 * t51 * t52 * t53 * t65 * 2.0)
				+ gz
						* (t79 + t90 + t91
								- b_wy * t40 * t49 * t51 * t52 * t53 * 2.0
								- b_wx * b_wz * t40 * t51 * t62 * t63 * 4.0
								- b_wy * dt * t40 * t51 * t54 * t62)
				+ gy
						* (t71 + t72 + t73 + t74
								- b_wx * b_wy * t40 * t51 * t62 * t63 * 4.0
								- b_wz * dt * t40 * t51 * t54 * t66);
		F(17,12) = gz
				* (t67 + t95 + t96 - b_wy * t41 * t49 * t52 * t53 * t55 * 2.0
						- b_wx * b_wz * t41 * t55 * t62 * t63 * 4.0
						- b_wy * dt * t41 * t54 * t55 * t62)
				- gx
						* (t75 + t92 + t94 - t41 * t55 * t62 * t63 * t64 * 4.0
								- t41 * t55 * t62 * t63 * t65 * 4.0)
				+ gy
						* (t68 + t76 + t77 + t78
								- b_wx * b_wy * t41 * t55 * t62 * t63 * 4.0
								- b_wz * dt * t41 * t54 * t55 * t66);
		F(17,13) = -gy
				* (t67 + t82 + t83 - b_wz * t42 * t49 * t52 * t53 * t56 * 2.0
						- b_wz * dt * t42 * t54 * t56 * t62
						- b_wx * b_wy * dt * t42 * t49 * t52 * t53 * t56 * 2.0)
				+ gz
						* (t68 + t102 + t103
								- b_wy * t42 * t49 * t52 * t53 * t56 * 2.0
								- b_wx * b_wz * t42 * t56 * t62 * t63 * 4.0
								- b_wy * dt * t42 * t54 * t56 * t62)
				- gx
						* (t80 + t81 + t101 - t42 * t56 * t62 * t63 * t64 * 4.0
								- t42 * t56 * t62 * t63 * t65 * 4.0);
		F(17,17) = t54 * t62 * t64 * -2.0 - t54 * t62 * t65 * 2.0 + 1.0;
		F(17,18) = t84 - b_wz * t49 * t50 * t53 * 2.0;
		F(17,19) = t107 + t108;
		F(18,11) = -gy
				* (-t69 + t85 + t87 - t40 * t51 * t62 * t63 * t70 * 4.0
						+ dt * t40 * t49 * t51 * t52 * t53 * t65 * 2.0)
				- gz
						* (t67 + t88 + t89
								- b_wx * t40 * t49 * t51 * t52 * t53 * 2.0
								- b_wx * dt * t40 * t51 * t54 * t62
								- b_wy * b_wz * dt * t40 * t49 * t51 * t52 * t53
										* 2.0)
				+ gx
						* (t71 - t72 - t73 + t74
								- b_wx * b_wy * t40 * t51 * t62 * t63 * 4.0
								+ b_wz * dt * t40 * t51 * t54 * t66);
		F(18,12) = -gy
				* (t75 + t93 - t41 * t55 * t62 * t63 * t65 * 4.0
						- t41 * t55 * t62 * t63 * t70 * 4.0)
				+ gx
						* (t68 - t76 - t77 + t78
								- b_wx * b_wy * t41 * t55 * t62 * t63 * 4.0
								+ b_wz * dt * t41 * t54 * t55 * t66)
				+ gz
						* (t79 + t97 + t98 + t99
								- b_wy * b_wz * t41 * t55 * t62 * t63 * 4.0
								- b_wx * dt * t41 * t54 * t55 * t66);
		F(18,13) = gx
				* (t67 + t82 - t83 - b_wz * t42 * t49 * t52 * t53 * t56 * 2.0
						- b_wz * dt * t42 * t54 * t56 * t62
						+ b_wx * b_wy * dt * t42 * t49 * t52 * t53 * t56 * 2.0)
				- gy
						* (t80 + t81 + t100 - t42 * t56 * t62 * t63 * t65 * 4.0
								- t42 * t56 * t62 * t63 * t70 * 4.0)
				+ gz
						* (t71 + t104 + t105 + t106
								- b_wy * b_wz * t42 * t56 * t62 * t63 * 4.0
								- b_wx * dt * t42 * t54 * t56 * t66);
		F(18,17) = t84 + b_wz * t49 * t50 * t53 * 2.0;
		F(18,18) = t54 * t62 * t65 * -2.0 - t54 * t62 * t70 * 2.0 + 1.0;
		F(18,19) = t109 - b_wx * t49 * t50 * t53 * 2.0;
		F(19,11) = -gz
				* (t85 - t86 + t87 - t40 * t51 * t62 * t63 * t70 * 4.0
						+ dt * t40 * t49 * t51 * t52 * t53 * t64 * 2.0)
				+ gy
						* (t67 + t88 - t89
								- b_wx * t40 * t49 * t51 * t52 * t53 * 2.0
								- b_wx * dt * t40 * t51 * t54 * t62
								+ b_wy * b_wz * dt * t40 * t49 * t51 * t52 * t53
										* 2.0)
				+ gx
						* (t79 - t90 + t91
								+ b_wy * t40 * t49 * t51 * t52 * t53 * 2.0
								- b_wx * b_wz * t40 * t51 * t62 * t63 * 4.0
								+ b_wy * dt * t40 * t51 * t54 * t62);
		F(19,12) = gy
				* (t79 - t97 - t98 + t99
						- b_wy * b_wz * t41 * t55 * t62 * t63 * 4.0
						+ b_wx * dt * t41 * t54 * t55 * t66)
				- gz
						* (t92 + t93 + t94 - t41 * t55 * t62 * t63 * t64 * 4.0
								- t41 * t55 * t62 * t63 * t70 * 4.0)
				- gx
						* (t67 + t95 - t96
								- b_wy * t41 * t49 * t52 * t53 * t55 * 2.0
								+ b_wx * b_wz * t41 * t55 * t62 * t63 * 4.0
								- b_wy * dt * t41 * t54 * t55 * t62);
		F(19,13) = -gz
				* (t100 + t101 - t42 * t56 * t62 * t63 * t64 * 4.0
						- t42 * t56 * t62 * t63 * t70 * 4.0)
				+ gy
						* (t71 - t104 - t105 + t106
								- b_wy * b_wz * t42 * t56 * t62 * t63 * 4.0
								+ b_wx * dt * t42 * t54 * t56 * t66)
				+ gx
						* (t68 - t102 + t103
								+ b_wy * t42 * t49 * t52 * t53 * t56 * 2.0
								- b_wx * b_wz * t42 * t56 * t62 * t63 * 4.0
								+ b_wy * dt * t42 * t54 * t56 * t62);
		F(19,17) = t107 - t108;
		F(19,18) = t109 + b_wx * t49 * t50 * t53 * 2.0;
		F(19,19) = t54 * t62 * t64 * -2.0 - t54 * t62 * t70 * 2.0 + 1.0;
		F(20,20) = 1.0;
		F(21,21) = 1.0;
		F(22,22) = 1.0;
		F(23,23) = 1.0;
		F(24,24) = 1.0;
		F(25,25) = 1.0;
	} else {
		F(0,0) = 1.0;
		F(0,3) = t34 - lambda * qy * t7 * 2.0;
		F(0,4) = t23 + t36;
		F(0,5) = t31 - lambda * qw * t7 * 2.0 - lambda * qy * t10 * 4.0;
		F(0,6) = t22 + t33 - lambda * qz * t10 * 4.0;
		F(0,7) = -t7 * t17 + t5 * t20 - t10 * t15;
		F(0,8) = -dt * lambda * t15;
		F(0,9) = -dt * lambda * t17 + dt * lambda * t20;
		F(0,14) = lambda * t3 * t15 * (-1.0 / 2.0);
		F(0,15) = lambda * t3 * t20 * (1.0 / 2.0);
		F(0,16) = lambda * t3 * t17 * (-1.0 / 2.0);
		F(1,1) = 1.0;
		F(1,3) = t22 - lambda * qz * t10 * 2.0;
		F(1,4) = t32 + lambda * qw * t7 * 2.0 - lambda * qx * t5 * 4.0;
		F(1,5) = t23 + t35;
		F(1,6) = lambda * qw * t10 * -2.0 + lambda * qy * t7 * 2.0
				- lambda * qz * t5 * 4.0;
		F(1,7) = t7 * t27 - t10 * t24 - t5 * t30;
		F(1,8) = -dt * lambda * t24;
		F(1,9) = dt * lambda * t27 - dt * lambda * t30;
		F(1,14) = lambda * t3 * t24 * (-1.0 / 2.0);
		F(1,15) = lambda * t3 * t30 * (-1.0 / 2.0);
		F(1,16) = lambda * t3 * t27 * (1.0 / 2.0);
		F(2,2) = 1.0;
		F(2,3) = -t31 + t32;
		F(2,4) = -t33 - lambda * qx * t7 * 4.0 + lambda * qz * t10 * 2.0;
		F(2,5) = t34 + lambda * qw * t10 * 2.0 - lambda * qy * t7 * 4.0;
		F(2,6) = t35 + t36;
		F(2,7) = -t5 * t38 - t7 * t39 + t10 * t37;
		F(2,8) = dt * lambda * t37;
		F(2,9) = -dt * lambda * t38 - dt * lambda * t39;
		F(2,14) = lambda * t3 * t37 * (1.0 / 2.0);
		F(2,15) = lambda * t3 * t38 * (-1.0 / 2.0);
		F(2,16) = lambda * t3 * t39 * (-1.0 / 2.0);
		F(3,3) = 1.0;
		F(4,4) = 1.0;
		F(5,5) = 1.0;
		F(6,6) = 1.0;
		F(7,7) = 1.0;
		F(8,8) = 1.0;
		F(8,14) = dt;
		F(9,9) = 1.0;
		F(9,15) = dt;
		F(10,9) = 1.0;
		F(10,16) = dt;
		F(11,11) = 1.0;
		F(12,12) = 1.0;
		F(13,13) = 1.0;
		F(14,14) = 1.0;
		F(15,15) = 1.0;
		F(16,16) = 1.0;
		F(17,17) = 1.0;
		F(18,18) = 1.0;
		F(19,19) = 1.0;
		F(20,20) = 1.0;
		F(21,21) = 1.0;
		F(22,22) = 1.0;
		F(23,23) = 1.0;
		F(24,24) = 1.0;
		F(25,25) = 1.0;
	}

}

