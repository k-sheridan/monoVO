/*
 * PoseEKF.cpp
 *
 *  Created on: Oct 31, 2017
 *      Author: kevin
 */

#include <PoseEKF.h>

PoseEKF::PoseEKF()
{
	// TODO Auto-generated constructor stub
	this->resetState();
}

PoseEKF::PoseEKF(ros::Time start)
{
	this->state.t = start;
	this->resetState();
}

PoseEKF::~PoseEKF()
{
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
	if (omega_norm > 1e-10)
	{ // this ensures that we don't divide by zero. im omega is zero that means that the delta quaternion is the identity quaternion
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
		Eigen::Matrix<double, 6, 6> cov, ros::Time t_measured)
{

	//we measure an euler angle, but our state stores its orientation in a quaternion
	Eigen::Vector3d omega = pose.log().template tail<3>();

	Eigen::Matrix<double, 4, 3> J;
	J.setZero();
	this->computeAngle2QuaternionJacobian(omega, J);

}

void PoseEKF::computeAngle2QuaternionJacobian(Eigen::Vector3d angle,
		Eigen::Matrix<double, 4, 3>& J)
{
	if (angle.squaredNorm() > 0.00000001)
	{
		  double t2 = fabs(angle.x());
		  double t3 = fabs(angle.y());
		  double t4 = fabs(angle.z());
		  double t5 = t2*t2;
		  double t6 = t3*t3;
		  double t7 = t4*t4;
		  double  t8 = t5+t6+t7;
		  double  t9 = sqrt(t8);
		  double  t10 = t9*(1.0/2.0);
		  double  t11 = sin(t10);
		  double  t12 = 1.0/sqrt(t8);
		  double  t13 = (angle.x()/fabs(angle.x()));
		  double  t14 = (angle.y()/fabs(angle.y()));
		  double   t15 = 1.0/pow(t8,3.0/2.0);
		  double  t16 = cos(t10);
		  double  t17 = 1.0/t8;
		  double  t18 = (angle.z()/fabs(angle.z()));
		  double  t19 = t11*t12;
		  J(0,0) = t2*t11*t12*t13*(-1.0/2.0);
		  J(0,1) = t3*t11*t12*t14*(-1.0/2.0);
		  J(0,2) = t4*t11*t12*t18*(-1.0/2.0);
		  J(1,0) = t19-angle.x()*t2*t11*t13*t15+angle.x()*t2*t13*t16*t17*(1.0/2.0);
		  J(1,1) = -angle.x()*t3*t11*t14*t15+angle.x()*t3*t14*t16*t17*(1.0/2.0);
		  J(1,2) = -angle.x()*t4*t11*t15*t18+angle.x()*t4*t16*t17*t18*(1.0/2.0);
		  J(2,0) = -angle.y()*t2*t11*t13*t15+angle.y()*t2*t13*t16*t17*(1.0/2.0);
		  J(2,1) = t19-angle.y()*t3*t11*t14*t15+angle.y()*t3*t14*t16*t17*(1.0/2.0);
		  J(2,2) = -angle.y()*t4*t11*t15*t18+angle.y()*t4*t16*t17*t18*(1.0/2.0);
		  J(3,0) = -angle.z()*t2*t11*t13*t15+angle.z()*t2*t13*t16*t17*(1.0/2.0);
		  J(3,1) = -angle.z()*t3*t11*t14*t15+angle.z()*t3*t14*t16*t17*(1.0/2.0);
		  J(3,2) = t19-angle.z()*t4*t11*t15*t18+angle.z()*t4*t16*t17*t18*(1.0/2.0);

	}
	else
	{
		double t2 = 1;
		double t3 = fabs(angle.x());
		double t4 = (angle.x()/fabs(angle.x()));
		double t5 = fabs(angle.y());
		double  t6 = fabs(angle.z());
		double  t7 = t2*t2;
		double  t8 = (angle.y()/fabs(angle.y()));
		double  t9 = t3*t3;
		double  t10 = t5*t5;
		double   t11 = t6*t6;
		double  t12 = t9+t10+t11;
		double  t13 = (angle.z()/fabs(angle.z()));
		double  t14 = t2*t3*t4*(1.0/2.4E1);
		double  t22 = t3*t4*t7*t12*(1.0/9.6E2);
		double  t15 = t14-t22;
		double   t16 = t12*t12;
		double   t17 = t7*t16*2.604166666666667E-4;
		double   t18 = t2*t5*t8*(1.0/2.4E1);
		double   t23 = t5*t7*t8*t12*(1.0/9.6E2);
		double   t19 = t18-t23;
		double   t20 = t2*t6*t13*(1.0/2.4E1);
		double   t24 = t6*t7*t12*t13*(1.0/9.6E2);
		double   t21 = t20-t24;
		  J(0,0) = -t2*t3*t4+t3*t4*t7*t12*(1.0/9.6E1);
		  J(0,1) = -t2*t5*t8+t5*t7*t8*t12*(1.0/9.6E1);
		  J(0,2) = -t2*t6*t13+t6*t7*t12*t13*(1.0/9.6E1);
		  J(1,0) = t17-angle.x()*t15-t2*t12*(1.0/4.8E1)+1.0/2.0;
		  J(1,1) = -angle.x()*t19;
		  J(1,2) = -angle.x()*t21;
		  J(2,0) = -angle.y()*t15;
		  J(2,1) = t17-angle.y()*t19-t2*t12*(1.0/4.8E1)+1.0/2.0;
		  J(2,2) = -angle.y()*t21;
		  J(3,0) = -angle.z()*t15;
		  J(3,1) = -angle.z()*t19;
		  J(3,2) = t17-angle.z()*t21-t2*t12*(1.0/4.8E1)+1.0/2.0;

	}
}

/*
 * due to the scale factor, and quaternion, and body frame velocity... this is not a linear process
 */
void PoseEKF::computeStateTransitionJacobian(State& from_state, double dt,
		Eigen::Matrix<double, STATE_SIZE, STATE_SIZE>& F)
{

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

	if (b_wx * b_wx + b_wy * b_wy + b_wz * b_wz > 0.00000001)
	{
		double t2 = b_dy * dt;
		double t3 = dt * dt;
		double t4 = b_ay * t3 * (1.0 / 2.0);
		double t5 = t2 + t4;
		double t6 = b_az * t3 * (1.0 / 2.0);
		double t7 = t2 + t6;
		double t8 = b_dx * dt;
		double t9 = b_ax * t3 * (1.0 / 2.0);
		double t10 = t8 + t9;
		double t11 = b_dy * 2.0;
		double t12 = qy * qy;
		double t13 = t12 * 2.0;
		double t14 = qz * qz;
		double t15 = t14 * 2.0;
		double t16 = t13 + t15 - 1.0;
		double t17 = qw * qy * 2.0;
		double t46 = qx * qz * 2.0;
		double t18 = t17 - t46;
		double t19 = qw * qz * 2.0;
		double t20 = qx * qy * 2.0;
		double t21 = t19 + t20;
		double t22 = lambda * qx * t7 * 2.0;
		double t23 = lambda * qz * t7 * 2.0;
		double t24 = b_ay * dt;
		double t25 = t11 + t24;
		double t26 = b_dx * 2.0;
		double t27 = b_ax * dt;
		double t28 = t26 + t27;
		double t29 = b_az * dt;
		double t30 = t11 + t29;
		double t31 = qw * qz;
		double t32 = qx * qy;
		double t33 = qw * qx * 2.0;
		double t34 = qy * qz * 2.0;
		double t35 = t33 + t34;
		double t36 = qx * qx;
		double t37 = t36 * 2.0;
		double t38 = t15 + t37 - 1.0;
		double t39 = t31 - t32;
		double t40 = lambda * qx * t5 * 2.0;
		double t41 = lambda * qy * t10 * 2.0;
		double t42 = lambda * qw * t5 * 2.0;
		double t43 = lambda * qz * t5 * 2.0;
		double t44 = lambda * qx * t10 * 2.0;
		double t45 = lambda * qy * t5 * 2.0;
		double t47 = t17 + t46;
		double t48 = t33 - t34;
		double t49 = t13 + t37 - 1.0;
		double t50 = qw * qy;
		double t51 = qw * qx;
		double t52 = qy * qz;
		double t53 = fabs(b_wx);
		double t54 = fabs(b_wy);
		double t55 = fabs(b_wz);
		double t56 = t53 * t53;
		double t57 = t54 * t54;
		double t58 = t55 * t55;
		double t59 = t56 + t57 + t58;
		double t60 = sqrt(t59);
		double t61 = dt * t60 * (1.0 / 2.0);
		double t62 = sin(t61);
		double t63 = 1.0 / sqrt(t59);
		double t64 = (b_wx / fabs(b_wx));
		double t65 = 1.0 / pow(t59, 3.0 / 2.0);
		double t66 = cos(t61);
		double t67 = 1.0 / t59;
		double t68 = (b_wy / fabs(b_wy));
		double t69 = (b_wz / fabs(b_wz));
		double t70 = b_wz * t62 * t63;
		double t71 = b_wx * t62 * t63;
		double t72 = qw * t62 * t63;
		double t73 = b_wy * t62 * t63;
		double t74 = qy * t62 * t63;
		double t75 = t62 * t62;
		double t76 = 1.0 / (t59 * t59);
		double t77 = b_wy * b_wy;
		double t78 = b_wz * b_wz;
		double t79 = t66 * t66;
		double t80 = dt * t60;
		double t81 = sin(t80);
		double t82 = t63 * t81;
		double t83 = b_wx * t67 * t75 * 2.0;
		double t84 = t53 * t64 * t75 * t76 * t78 * 4.0;
		double t85 = b_wx * b_wx;
		double t86 = b_wy * t67 * t75 * 2.0;
		double t87 = b_wz * t53 * t62 * t64 * t65 * t66 * 2.0;
		double t88 = b_wz * dt * t53 * t64 * t67 * t75;
		double t89 = b_wx * b_wy * dt * t53 * t62 * t64 * t65 * t66 * 2.0;
		double t90 = dt * t54 * t62 * t65 * t66 * t68 * t78 * 2.0;
		double t91 = b_wz * t54 * t62 * t65 * t66 * t68 * 2.0;
		double t92 = b_wz * dt * t54 * t67 * t68 * t75;
		double t93 = b_wx * b_wy * dt * t54 * t62 * t65 * t66 * t68 * 2.0;
		double t94 = b_wz * t67 * t75 * 2.0;
		double t95 = b_wz * t67 * t75 * 4.0;
		double t96 = dt * t55 * t62 * t65 * t66 * t69 * t78 * 2.0;
		double t97 = b_wz * dt * t55 * t67 * t69 * t79;
		double t98 = b_wx * b_wy * t55 * t69 * t75 * t76 * 4.0;
		double t99 = b_wx * b_wy * t67 * t75 * 2.0;
		double t100 = b_wx * t67 * t75 * 4.0;
		double t101 = t53 * t64 * t75 * t76 * t77 * 4.0;
		double t102 = dt * t53 * t62 * t64 * t65 * t66 * t85 * 2.0;
		double t103 = b_wx * dt * t53 * t64 * t67 * t79;
		double t104 = b_wy * b_wz * t53 * t64 * t75 * t76 * 4.0;
		double t105 = b_wy * dt * t53 * t64 * t67 * t79;
		double t106 = b_wx * b_wz * dt * t53 * t62 * t64 * t65 * t66 * 2.0;
		double t107 = b_wy * t67 * t75 * 4.0;
		double t108 = dt * t54 * t62 * t65 * t66 * t68 * t85 * 2.0;
		double t109 = dt * t54 * t62 * t65 * t66 * t68 * t77 * 2.0;
		double t110 = b_wy * dt * t54 * t67 * t68 * t79;
		double t111 = b_wx * b_wz * dt * t54 * t62 * t65 * t66 * t68 * 2.0;
		double t112 = b_wx * t54 * t62 * t65 * t66 * t68 * 2.0;
		double t113 = b_wx * dt * t54 * t67 * t68 * t75;
		double t114 = b_wy * b_wz * dt * t54 * t62 * t65 * t66 * t68 * 2.0;
		double t115 = dt * t55 * t62 * t65 * t66 * t69 * t85 * 2.0;
		double t116 = dt * t55 * t62 * t65 * t66 * t69 * t77 * 2.0;
		double t117 = b_wy * dt * t55 * t67 * t69 * t79;
		double t118 = b_wx * b_wz * dt * t55 * t62 * t65 * t66 * t69 * 2.0;
		double t119 = b_wx * t55 * t62 * t65 * t66 * t69 * 2.0;
		double t120 = b_wx * dt * t55 * t67 * t69 * t75;
		double t121 = b_wy * b_wz * dt * t55 * t62 * t65 * t66 * t69 * 2.0;
		double t122 = b_wy * t63 * t81;
		double t123 = b_wx * b_wz * t67 * t75 * 2.0;
		double t124 = b_wy * b_wz * t67 * t75 * 2.0;
		F(0, 0) = 1.0;
		F(0, 3) = t43 - lambda * qy * t7 * 2.0;
		F(0, 4) = t23 + t45;
		F(0, 5) = t40 - lambda * qw * t7 * 2.0 - lambda * qy * t10 * 4.0;
		F(0, 6) = t22 + t42 - lambda * qz * t10 * 4.0;
		F(0, 7) = dt * t16 * t28 * (-1.0 / 2.0) + dt * t21 * t25 * (1.0 / 2.0)
				- dt * t18 * t30 * (1.0 / 2.0);
		F(0, 8) = -dt * lambda * t16;
		F(0, 9) = -dt * lambda * t18 + dt * lambda * t21;
		F(0, 14) = lambda * t3 * t16 * (-1.0 / 2.0);
		F(0, 15) = lambda * t3 * (t31 + t32);
		F(0, 16) = -lambda * t3 * (t50 - qx * qz);
		F(1, 1) = 1.0;
		F(1, 3) = t22 - lambda * qz * t10 * 2.0;
		F(1, 4) = t41 + lambda * qw * t7 * 2.0 - lambda * qx * t5 * 4.0;
		F(1, 5) = t23 + t44;
		F(1, 6) = lambda * qw * t10 * -2.0 + lambda * qy * t7 * 2.0
				- lambda * qz * t5 * 4.0;
		F(1, 7) = dt * t25 * t38 * (-1.0 / 2.0) + dt * t30 * t35 * (1.0 / 2.0)
				- dt * t28 * (t19 - t20) * (1.0 / 2.0);
		F(1, 8) = dt * lambda * t39 * -2.0;
		F(1, 9) = dt * lambda * t35 - dt * lambda * t38;
		F(1, 14) = -lambda * t3 * t39;
		F(1, 15) = lambda * t3 * t38 * (-1.0 / 2.0);
		F(1, 16) = lambda * t3 * (t51 + t52);
		F(2, 2) = 1.0;
		F(2, 3) = -t40 + t41;
		F(2, 4) = -t42 - lambda * qx * t7 * 4.0 + lambda * qz * t10 * 2.0;
		F(2, 5) = t43 + lambda * qw * t10 * 2.0 - lambda * qy * t7 * 4.0;
		F(2, 6) = t44 + t45;
		F(2, 7) = dt * t25 * t48 * (-1.0 / 2.0) + dt * t28 * t47 * (1.0 / 2.0)
				- dt * t30 * t49 * (1.0 / 2.0);
		F(2, 8) = dt * lambda * t47;
		F(2, 9) = -dt * lambda * t48 - dt * lambda * t49;
		F(2, 14) = lambda * t3 * (t50 + qx * qz);
		F(2, 15) = -lambda * t3 * (t51 - t52);
		F(2, 16) = lambda * t3 * t49 * (-1.0 / 2.0);
		F(3, 3) = t66;
		F(3, 4) = -b_wx * t62 * t63;
		F(3, 5) = -b_wy * t62 * t63;
		F(3, 6) = -b_wz * t62 * t63;
		F(3, 11) = -qx * t62 * t63 + b_wx * qx * t53 * t62 * t64 * t65
				+ b_wy * qy * t53 * t62 * t64 * t65
				+ b_wz * qz * t53 * t62 * t64 * t65
				- dt * qw * t53 * t62 * t63 * t64 * (1.0 / 2.0)
				- b_wx * dt * qx * t53 * t64 * t66 * t67 * (1.0 / 2.0)
				- b_wy * dt * qy * t53 * t64 * t66 * t67 * (1.0 / 2.0)
				- b_wz * dt * qz * t53 * t64 * t66 * t67 * (1.0 / 2.0);
		F(3, 12) = -qy * t62 * t63 + b_wx * qx * t54 * t62 * t65 * t68
				+ b_wy * qy * t54 * t62 * t65 * t68
				+ b_wz * qz * t54 * t62 * t65 * t68
				- dt * qw * t54 * t62 * t63 * t68 * (1.0 / 2.0)
				- b_wx * dt * qx * t54 * t66 * t67 * t68 * (1.0 / 2.0)
				- b_wy * dt * qy * t54 * t66 * t67 * t68 * (1.0 / 2.0)
				- b_wz * dt * qz * t54 * t66 * t67 * t68 * (1.0 / 2.0);
		F(3, 13) = -qz * t62 * t63 + b_wx * qx * t55 * t62 * t65 * t69
				+ b_wy * qy * t55 * t62 * t65 * t69
				+ b_wz * qz * t55 * t62 * t65 * t69
				- dt * qw * t55 * t62 * t63 * t69 * (1.0 / 2.0)
				- b_wx * dt * qx * t55 * t66 * t67 * t69 * (1.0 / 2.0)
				- b_wy * dt * qy * t55 * t66 * t67 * t69 * (1.0 / 2.0)
				- b_wz * dt * qz * t55 * t66 * t67 * t69 * (1.0 / 2.0);
		F(4, 3) = t71;
		F(4, 4) = t66;
		F(4, 5) = t70;
		F(4, 6) = -b_wy * t62 * t63;
		F(4, 11) = t72 - b_wx * qw * t53 * t62 * t64 * t65
				- b_wz * qy * t53 * t62 * t64 * t65
				+ b_wy * qz * t53 * t62 * t64 * t65
				- dt * qx * t53 * t62 * t63 * t64 * (1.0 / 2.0)
				+ b_wx * dt * qw * t53 * t64 * t66 * t67 * (1.0 / 2.0)
				+ b_wz * dt * qy * t53 * t64 * t66 * t67 * (1.0 / 2.0)
				- b_wy * dt * qz * t53 * t64 * t66 * t67 * (1.0 / 2.0);
		F(4, 12) = -qz * t62 * t63 - b_wx * qw * t54 * t62 * t65 * t68
				- b_wz * qy * t54 * t62 * t65 * t68
				+ b_wy * qz * t54 * t62 * t65 * t68
				- dt * qx * t54 * t62 * t63 * t68 * (1.0 / 2.0)
				+ b_wx * dt * qw * t54 * t66 * t67 * t68 * (1.0 / 2.0)
				+ b_wz * dt * qy * t54 * t66 * t67 * t68 * (1.0 / 2.0)
				- b_wy * dt * qz * t54 * t66 * t67 * t68 * (1.0 / 2.0);
		F(4, 13) = t74 - b_wx * qw * t55 * t62 * t65 * t69
				- b_wz * qy * t55 * t62 * t65 * t69
				+ b_wy * qz * t55 * t62 * t65 * t69
				- dt * qx * t55 * t62 * t63 * t69 * (1.0 / 2.0)
				+ b_wx * dt * qw * t55 * t66 * t67 * t69 * (1.0 / 2.0)
				+ b_wz * dt * qy * t55 * t66 * t67 * t69 * (1.0 / 2.0)
				- b_wy * dt * qz * t55 * t66 * t67 * t69 * (1.0 / 2.0);
		F(5, 3) = t73;
		F(5, 4) = -t70;
		F(5, 5) = t66;
		F(5, 6) = t71;
		F(5, 11) = qz * t62 * t63 - b_wy * qw * t53 * t62 * t64 * t65
				+ b_wz * qx * t53 * t62 * t64 * t65
				- b_wx * qz * t53 * t62 * t64 * t65
				- dt * qy * t53 * t62 * t63 * t64 * (1.0 / 2.0)
				+ b_wy * dt * qw * t53 * t64 * t66 * t67 * (1.0 / 2.0)
				- b_wz * dt * qx * t53 * t64 * t66 * t67 * (1.0 / 2.0)
				+ b_wx * dt * qz * t53 * t64 * t66 * t67 * (1.0 / 2.0);
		F(5, 12) = t72 - b_wy * qw * t54 * t62 * t65 * t68
				+ b_wz * qx * t54 * t62 * t65 * t68
				- b_wx * qz * t54 * t62 * t65 * t68
				- dt * qy * t54 * t62 * t63 * t68 * (1.0 / 2.0)
				+ b_wy * dt * qw * t54 * t66 * t67 * t68 * (1.0 / 2.0)
				- b_wz * dt * qx * t54 * t66 * t67 * t68 * (1.0 / 2.0)
				+ b_wx * dt * qz * t54 * t66 * t67 * t68 * (1.0 / 2.0);
		F(5, 13) = -qx * t62 * t63 - b_wy * qw * t55 * t62 * t65 * t69
				+ b_wz * qx * t55 * t62 * t65 * t69
				- b_wx * qz * t55 * t62 * t65 * t69
				- dt * qy * t55 * t62 * t63 * t69 * (1.0 / 2.0)
				+ b_wy * dt * qw * t55 * t66 * t67 * t69 * (1.0 / 2.0)
				- b_wz * dt * qx * t55 * t66 * t67 * t69 * (1.0 / 2.0)
				+ b_wx * dt * qz * t55 * t66 * t67 * t69 * (1.0 / 2.0);
		F(6, 3) = t70;
		F(6, 4) = t73;
		F(6, 5) = -t71;
		F(6, 6) = t66;
		F(6, 11) = -t74 - b_wz * qw * t53 * t62 * t64 * t65
				- b_wy * qx * t53 * t62 * t64 * t65
				+ b_wx * qy * t53 * t62 * t64 * t65
				- dt * qz * t53 * t62 * t63 * t64 * (1.0 / 2.0)
				+ b_wz * dt * qw * t53 * t64 * t66 * t67 * (1.0 / 2.0)
				+ b_wy * dt * qx * t53 * t64 * t66 * t67 * (1.0 / 2.0)
				- b_wx * dt * qy * t53 * t64 * t66 * t67 * (1.0 / 2.0);
		F(6, 12) = qx * t62 * t63 - b_wz * qw * t54 * t62 * t65 * t68
				- b_wy * qx * t54 * t62 * t65 * t68
				+ b_wx * qy * t54 * t62 * t65 * t68
				- dt * qz * t54 * t62 * t63 * t68 * (1.0 / 2.0)
				+ b_wz * dt * qw * t54 * t66 * t67 * t68 * (1.0 / 2.0)
				+ b_wy * dt * qx * t54 * t66 * t67 * t68 * (1.0 / 2.0)
				- b_wx * dt * qy * t54 * t66 * t67 * t68 * (1.0 / 2.0);
		F(6, 13) = t72 - b_wz * qw * t55 * t62 * t65 * t69
				- b_wy * qx * t55 * t62 * t65 * t69
				+ b_wx * qy * t55 * t62 * t65 * t69
				- dt * qz * t55 * t62 * t63 * t69 * (1.0 / 2.0)
				+ b_wz * dt * qw * t55 * t66 * t67 * t69 * (1.0 / 2.0)
				+ b_wy * dt * qx * t55 * t66 * t67 * t69 * (1.0 / 2.0)
				- b_wx * dt * qy * t55 * t66 * t67 * t69 * (1.0 / 2.0);
		F(7, 7) = 1.0;
		F(8, 8) = 1.0;
		F(8, 14) = dt;
		F(9, 9) = 1.0;
		F(9, 15) = dt;
		F(10, 9) = 1.0;
		F(10, 16) = dt;
		F(11, 11) = 1.0;
		F(12, 12) = 1.0;
		F(13, 13) = 1.0;
		F(14, 14) = 1.0;
		F(15, 15) = 1.0;
		F(16, 16) = 1.0;
		F(17, 11) = gx
				* (t84 + t101 - dt * t53 * t62 * t64 * t65 * t66 * t77 * 2.0
						- dt * t53 * t62 * t64 * t65 * t66 * t78 * 2.0)
				+ gz
						* (t94 + t105 + t106
								- b_wy * t53 * t62 * t64 * t65 * t66 * 2.0
								- b_wx * b_wz * t53 * t64 * t75 * t76 * 4.0
								- b_wy * dt * t53 * t64 * t67 * t75)
				+ gy
						* (t86 + t87 + t88 + t89
								- b_wx * b_wy * t53 * t64 * t75 * t76 * 4.0
								- b_wz * dt * t53 * t64 * t67 * t79);
		F(17, 12) = gz
				* (t82 + t110 + t111 - b_wy * t54 * t62 * t65 * t66 * t68 * 2.0
						- b_wx * b_wz * t54 * t68 * t75 * t76 * 4.0
						- b_wy * dt * t54 * t67 * t68 * t75)
				- gx
						* (t90 + t107 + t109 - t54 * t68 * t75 * t76 * t77 * 4.0
								- t54 * t68 * t75 * t76 * t78 * 4.0)
				+ gy
						* (t83 + t91 + t92 + t93
								- b_wx * b_wy * t54 * t68 * t75 * t76 * 4.0
								- b_wz * dt * t54 * t67 * t68 * t79);
		F(17, 13) = -gy
				* (t82 + t97 + t98 - b_wz * t55 * t62 * t65 * t66 * t69 * 2.0
						- b_wz * dt * t55 * t67 * t69 * t75
						- b_wx * b_wy * dt * t55 * t62 * t65 * t66 * t69 * 2.0)
				+ gz
						* (t83 + t117 + t118
								- b_wy * t55 * t62 * t65 * t66 * t69 * 2.0
								- b_wx * b_wz * t55 * t69 * t75 * t76 * 4.0
								- b_wy * dt * t55 * t67 * t69 * t75)
				- gx
						* (t95 + t96 + t116 - t55 * t69 * t75 * t76 * t77 * 4.0
								- t55 * t69 * t75 * t76 * t78 * 4.0);
		F(17, 17) = t67 * t75 * t77 * -2.0 - t67 * t75 * t78 * 2.0 + 1.0;
		F(17, 18) = t99 - b_wz * t63 * t81;
		F(17, 19) = t122 + t123;
		F(18, 11) = -gy
				* (-t84 + t100 + t102 - t53 * t64 * t75 * t76 * t85 * 4.0
						+ dt * t53 * t62 * t64 * t65 * t66 * t78 * 2.0)
				- gz
						* (t82 + t103 + t104
								- b_wx * t53 * t62 * t64 * t65 * t66 * 2.0
								- b_wx * dt * t53 * t64 * t67 * t75
								- b_wy * b_wz * dt * t53 * t62 * t64 * t65 * t66
										* 2.0)
				+ gx
						* (t86 - t87 - t88 + t89
								- b_wx * b_wy * t53 * t64 * t75 * t76 * 4.0
								+ b_wz * dt * t53 * t64 * t67 * t79);
		F(18, 12) = -gy
				* (t90 + t108 - t54 * t68 * t75 * t76 * t78 * 4.0
						- t54 * t68 * t75 * t76 * t85 * 4.0)
				+ gx
						* (t83 - t91 - t92 + t93
								- b_wx * b_wy * t54 * t68 * t75 * t76 * 4.0
								+ b_wz * dt * t54 * t67 * t68 * t79)
				+ gz
						* (t94 + t112 + t113 + t114
								- b_wy * b_wz * t54 * t68 * t75 * t76 * 4.0
								- b_wx * dt * t54 * t67 * t68 * t79);
		F(18, 13) = gx
				* (t82 + t97 - t98 - b_wz * t55 * t62 * t65 * t66 * t69 * 2.0
						- b_wz * dt * t55 * t67 * t69 * t75
						+ b_wx * b_wy * dt * t55 * t62 * t65 * t66 * t69 * 2.0)
				- gy
						* (t95 + t96 + t115 - t55 * t69 * t75 * t76 * t78 * 4.0
								- t55 * t69 * t75 * t76 * t85 * 4.0)
				+ gz
						* (t86 + t119 + t120 + t121
								- b_wy * b_wz * t55 * t69 * t75 * t76 * 4.0
								- b_wx * dt * t55 * t67 * t69 * t79);
		F(18, 17) = t99 + b_wz * t63 * t81;
		F(18, 18) = t67 * t75 * t78 * -2.0 - t67 * t75 * t85 * 2.0 + 1.0;
		F(18, 19) = t124 - b_wx * t63 * t81;
		F(19, 11) = -gz
				* (t100 - t101 + t102 - t53 * t64 * t75 * t76 * t85 * 4.0
						+ dt * t53 * t62 * t64 * t65 * t66 * t77 * 2.0)
				+ gy
						* (t82 + t103 - t104
								- b_wx * t53 * t62 * t64 * t65 * t66 * 2.0
								- b_wx * dt * t53 * t64 * t67 * t75
								+ b_wy * b_wz * dt * t53 * t62 * t64 * t65 * t66
										* 2.0)
				+ gx
						* (t94 - t105 + t106
								+ b_wy * t53 * t62 * t64 * t65 * t66 * 2.0
								- b_wx * b_wz * t53 * t64 * t75 * t76 * 4.0
								+ b_wy * dt * t53 * t64 * t67 * t75);
		F(19, 12) = gy
				* (t94 - t112 - t113 + t114
						- b_wy * b_wz * t54 * t68 * t75 * t76 * 4.0
						+ b_wx * dt * t54 * t67 * t68 * t79)
				- gz
						* (t107 + t108 + t109
								- t54 * t68 * t75 * t76 * t77 * 4.0
								- t54 * t68 * t75 * t76 * t85 * 4.0)
				- gx
						* (t82 + t110 - t111
								- b_wy * t54 * t62 * t65 * t66 * t68 * 2.0
								+ b_wx * b_wz * t54 * t68 * t75 * t76 * 4.0
								- b_wy * dt * t54 * t67 * t68 * t75);
		F(19, 13) = -gz
				* (t115 + t116 - t55 * t69 * t75 * t76 * t77 * 4.0
						- t55 * t69 * t75 * t76 * t85 * 4.0)
				+ gy
						* (t86 - t119 - t120 + t121
								- b_wy * b_wz * t55 * t69 * t75 * t76 * 4.0
								+ b_wx * dt * t55 * t67 * t69 * t79)
				+ gx
						* (t83 - t117 + t118
								+ b_wy * t55 * t62 * t65 * t66 * t69 * 2.0
								- b_wx * b_wz * t55 * t69 * t75 * t76 * 4.0
								+ b_wy * dt * t55 * t67 * t69 * t75);
		F(19, 17) = -t122 + t123;
		F(19, 18) = t124 + b_wx * t63 * t81;
		F(19, 19) = t67 * t75 * t77 * -2.0 - t67 * t75 * t85 * 2.0 + 1.0;
		F(20, 20) = 1.0;
		F(21, 21) = 1.0;
		F(22, 22) = 1.0;
		F(23, 23) = 1.0;
		F(24, 24) = 1.0;
		F(25, 25) = 1.0;

	}
	else
	{
		double t2 = b_dy * dt;
		double t3 = dt * dt;
		double t4 = b_ay * t3 * (1.0 / 2.0);
		double t5 = t2 + t4;
		double t6 = b_az * t3 * (1.0 / 2.0);
		double t7 = t2 + t6;
		double t8 = b_dx * dt;
		double t9 = b_ax * t3 * (1.0 / 2.0);
		double t10 = t8 + t9;
		double t11 = b_dy * 2.0;
		double t12 = qy * qy;
		double t13 = t12 * 2.0;
		double t14 = qz * qz;
		double t15 = t14 * 2.0;
		double t16 = t13 + t15 - 1.0;
		double t17 = qw * qy * 2.0;
		double t46 = qx * qz * 2.0;
		double t18 = t17 - t46;
		double t19 = qw * qz * 2.0;
		double t20 = qx * qy * 2.0;
		double t21 = t19 + t20;
		double t22 = lambda * qx * t7 * 2.0;
		double t23 = lambda * qz * t7 * 2.0;
		double t24 = b_ay * dt;
		double t25 = t11 + t24;
		double t26 = b_dx * 2.0;
		double t27 = b_ax * dt;
		double t28 = t26 + t27;
		double t29 = b_az * dt;
		double t30 = t11 + t29;
		double t31 = qw * qz;
		double t32 = qx * qy;
		double t33 = qw * qx * 2.0;
		double t34 = qy * qz * 2.0;
		double t35 = t33 + t34;
		double t36 = qx * qx;
		double t37 = t36 * 2.0;
		double t38 = t15 + t37 - 1.0;
		double t39 = t31 - t32;
		double t40 = lambda * qx * t5 * 2.0;
		double t41 = lambda * qy * t10 * 2.0;
		double t42 = lambda * qw * t5 * 2.0;
		double t43 = lambda * qz * t5 * 2.0;
		double t44 = lambda * qx * t10 * 2.0;
		double t45 = lambda * qy * t5 * 2.0;
		double t47 = t17 + t46;
		double t48 = t33 - t34;
		double t49 = t13 + t37 - 1.0;
		double t50 = qw * qy;
		double t51 = qw * qx;
		double t52 = qy * qz;
		double t53 = fabs(b_wx);
		double t54 = fabs(b_wy);
		double t55 = fabs(b_wz);
		double t57 = t53 * t53;
		double t58 = t54 * t54;
		double t59 = t55 * t55;
		double t56 = t57 + t58 + t59;
		double t60 = t3 * t3;
		double t61 = t56 * t56;
		double t62 = t60 * t61 * 2.604166666666667E-4;
		double t64 = t3 * t56 * (1.0 / 4.8E1);
		double t63 = t62 - t64 + 1.0 / 2.0;
		double t65 = (b_wx / fabs(b_wx));
		double t66 = t3 * t53 * t65 * (1.0 / 2.4E1);
		double t68 = t53 * t56 * t60 * t65 * (1.0 / 9.6E2);
		double t67 = t66 - t68;
		double t69 = (b_wy / fabs(b_wy));
		double t70 = t3 * t54 * t69 * (1.0 / 2.4E1);
		double t72 = t54 * t56 * t60 * t69 * (1.0 / 9.6E2);
		double t71 = t70 - t72;
		double t73 = (b_wz / fabs(b_wz));
		double t74 = t3 * t55 * t73 * (1.0 / 2.4E1);
		double t76 = t55 * t56 * t60 * t73 * (1.0 / 9.6E2);
		double t75 = t74 - t76;
		double t77 = t60 * t61 * (1.0 / 3.84E2);
		double t86 = t3 * t56 * (1.0 / 2.0);
		double t78 = t77 - t86 + 1.0;
		double t79 = t3 * t53 * t65;
		double t88 = t53 * t56 * t60 * t65 * (1.0 / 9.6E1);
		double t80 = t79 - t88;
		double t81 = t3 * t54 * t69;
		double t90 = t54 * t56 * t60 * t69 * (1.0 / 9.6E1);
		double t82 = t81 - t90;
		double t83 = t3 * t55 * t73;
		double t91 = t55 * t56 * t60 * t73 * (1.0 / 9.6E1);
		double t84 = t83 - t91;
		double t85 = b_wz * t63;
		double t87 = b_wx * t63;
		double t89 = qw * t63;
		double t92 = b_wy * t63;
		double t93 = qy * t63;
		double t94 = t63 * t63;
		double t95 = b_wy * b_wy;
		double t96 = b_wz * b_wz;
		double t97 = b_wx * t94 * 2.0;
		double t99 = t60 * t61;
		double t100 = t3 * t56 * 8.0E1;
		double t98 = t99 - t100 + 1.92E3;
		double t101 = t98 * t98;
		double t102 = t63 * t78 * 2.0;
		double t103 = b_wy * t94 * 2.0;
		double t104 = b_wz * t67 * t78 * 2.0;
		double t105 = b_wz * t63 * t80 * 2.0;
		double t106 = t63 * t67 * t96 * 4.0;
		double t107 = b_wz * t71 * t78 * 2.0;
		double t108 = b_wz * t63 * t82 * 2.0;
		double t109 = b_wx * b_wx;
		double t110 = t63 * t71 * t96 * 4.0;
		double t111 = b_wx * b_wy * t63 * t75 * 4.0;
		double t112 = t63 * t75 * t96 * 4.0;
		double t113 = b_wx * b_wy * t101 * 1.356336805555556E-7;
		double t114 = b_wy * b_wz * t63 * t67 * 4.0;
		double t115 = b_wz * t94 * 2.0;
		double t116 = b_wy * t67 * t78 * 2.0;
		double t117 = b_wy * t63 * t80 * 2.0;
		double t118 = b_wx * b_wz * t63 * t67 * 4.0;
		double t119 = t63 * t67 * t109 * 4.0;
		double t120 = t63 * t67 * t95 * 4.0;
		double t121 = b_wy * t71 * t78 * 2.0;
		double t122 = b_wy * t63 * t82 * 2.0;
		double t123 = b_wx * b_wz * t63 * t71 * 4.0;
		double t124 = b_wx * t71 * t78 * 2.0;
		double t125 = b_wx * t63 * t82 * 2.0;
		double t126 = t63 * t71 * t109 * 4.0;
		double t127 = t63 * t71 * t95 * 4.0;
		double t128 = b_wy * t75 * t78 * 2.0;
		double t129 = b_wy * t63 * t84 * 2.0;
		double t130 = b_wx * b_wz * t63 * t75 * 4.0;
		double t131 = b_wx * t75 * t78 * 2.0;
		double t132 = b_wx * t63 * t84 * 2.0;
		double t133 = t63 * t75 * t109 * 4.0;
		double t134 = t63 * t75 * t95 * 4.0;
		double t135 = b_wx * b_wz * t101 * 1.356336805555556E-7;
		double t136 = b_wy * t63 * t78 * 2.0;
		double t137 = b_wy * b_wz * t101 * 1.356336805555556E-7;
		F(0, 0) = 1.0;
		F(0, 3) = t43 - lambda * qy * t7 * 2.0;
		F(0, 4) = t23 + t45;
		F(0, 5) = t40 - lambda * qw * t7 * 2.0 - lambda * qy * t10 * 4.0;
		F(0, 6) = t22 + t42 - lambda * qz * t10 * 4.0;
		F(0, 7) = dt * t16 * t28 * (-1.0 / 2.0) + dt * t21 * t25 * (1.0 / 2.0)
				- dt * t18 * t30 * (1.0 / 2.0);
		F(0, 8) = -dt * lambda * t16;
		F(0, 9) = -dt * lambda * t18 + dt * lambda * t21;
		F(0, 14) = lambda * t3 * t16 * (-1.0 / 2.0);
		F(0, 15) = lambda * t3 * (t31 + t32);
		F(0, 16) = -lambda * t3 * (t50 - qx * qz);
		F(1, 1) = 1.0;
		F(1, 3) = t22 - lambda * qz * t10 * 2.0;
		F(1, 4) = t41 + lambda * qw * t7 * 2.0 - lambda * qx * t5 * 4.0;
		F(1, 5) = t23 + t44;
		F(1, 6) = lambda * qw * t10 * -2.0 + lambda * qy * t7 * 2.0
				- lambda * qz * t5 * 4.0;
		F(1, 7) = dt * t25 * t38 * (-1.0 / 2.0) + dt * t30 * t35 * (1.0 / 2.0)
				- dt * t28 * (t19 - t20) * (1.0 / 2.0);
		F(1, 8) = dt * lambda * t39 * -2.0;
		F(1, 9) = dt * lambda * t35 - dt * lambda * t38;
		F(1, 14) = -lambda * t3 * t39;
		F(1, 15) = lambda * t3 * t38 * (-1.0 / 2.0);
		F(1, 16) = lambda * t3 * (t51 + t52);
		F(2, 2) = 1.0;
		F(2, 3) = -t40 + t41;
		F(2, 4) = -t42 - lambda * qx * t7 * 4.0 + lambda * qz * t10 * 2.0;
		F(2, 5) = t43 + lambda * qw * t10 * 2.0 - lambda * qy * t7 * 4.0;
		F(2, 6) = t44 + t45;
		F(2, 7) = dt * t25 * t48 * (-1.0 / 2.0) + dt * t28 * t47 * (1.0 / 2.0)
				- dt * t30 * t49 * (1.0 / 2.0);
		F(2, 8) = dt * lambda * t47;
		F(2, 9) = -dt * lambda * t48 - dt * lambda * t49;
		F(2, 14) = lambda * t3 * (t50 + qx * qz);
		F(2, 15) = -lambda * t3 * (t51 - t52);
		F(2, 16) = lambda * t3 * t49 * (-1.0 / 2.0);
		F(3, 3) = t78;
		F(3, 4) = -b_wx * t63;
		F(3, 5) = -b_wy * t63;
		F(3, 6) = -b_wz * t63;
		F(3, 11) = -qw * t80 - qx * t63 + b_wx * qx * t67 + b_wy * qy * t67
				+ b_wz * qz * t67;
		F(3, 12) = -qw * t82 - qy * t63 + b_wx * qx * t71 + b_wy * qy * t71
				+ b_wz * qz * t71;
		F(3, 13) = -qw * t84 - qz * t63 + b_wx * qx * t75 + b_wy * qy * t75
				+ b_wz * qz * t75;
		F(4, 3) = t87;
		F(4, 4) = t78;
		F(4, 5) = t85;
		F(4, 6) = -b_wy * t63;
		F(4, 11) = t89 - qx * t80 - b_wx * qw * t67 - b_wz * qy * t67
				+ b_wy * qz * t67;
		F(4, 12) = -qx * t82 - qz * t63 - b_wx * qw * t71 - b_wz * qy * t71
				+ b_wy * qz * t71;
		F(4, 13) = t93 - qx * t84 - b_wx * qw * t75 - b_wz * qy * t75
				+ b_wy * qz * t75;
		F(5, 3) = t92;
		F(5, 4) = -t85;
		F(5, 5) = t78;
		F(5, 6) = t87;
		F(5, 11) = -qy * t80 + qz * t63 - b_wy * qw * t67 + b_wz * qx * t67
				- b_wx * qz * t67;
		F(5, 12) = t89 - qy * t82 - b_wy * qw * t71 + b_wz * qx * t71
				- b_wx * qz * t71;
		F(5, 13) = -qx * t63 - qy * t84 - b_wy * qw * t75 + b_wz * qx * t75
				- b_wx * qz * t75;
		F(6, 3) = t85;
		F(6, 4) = t92;
		F(6, 5) = -t87;
		F(6, 6) = t78;
		F(6, 11) = -t93 - qz * t80 - b_wz * qw * t67 - b_wy * qx * t67
				+ b_wx * qy * t67;
		F(6, 12) = qx * t63 - qz * t82 - b_wz * qw * t71 - b_wy * qx * t71
				+ b_wx * qy * t71;
		F(6, 13) = t89 - qz * t84 - b_wz * qw * t75 - b_wy * qx * t75
				+ b_wx * qy * t75;
		F(7, 7) = 1.0;
		F(8, 8) = 1.0;
		F(8, 14) = dt;
		F(9, 9) = 1.0;
		F(9, 15) = dt;
		F(10, 9) = 1.0;
		F(10, 16) = dt;
		F(11, 11) = 1.0;
		F(12, 12) = 1.0;
		F(13, 13) = 1.0;
		F(14, 14) = 1.0;
		F(15, 15) = 1.0;
		F(16, 16) = 1.0;
		F(17, 11) = gx * (t106 + t120)
				+ gy * (t103 + t104 + t105 - b_wx * b_wy * t63 * t67 * 4.0)
				- gz * (t116 + t117 + t118 - b_wz * t94 * 2.0);
		F(17, 12) = gx * (t110 + t127 - b_wy * t94 * 4.0)
				+ gy * (t97 + t107 + t108 - b_wx * b_wy * t63 * t71 * 4.0)
				- gz * (t121 + t122 + t123 - t63 * t78 * 2.0);
		F(17, 13) = gx * (t112 + t134 - b_wz * t94 * 4.0)
				- gz * (-t97 + t128 + t129 + t130)
				- gy
						* (t102 + t111 - b_wz * t63 * t84 * 2.0
								- b_wz * t75 * t78 * 2.0);
		F(17, 17) = t95 * t101 * (-1.356336805555556E-7)
				- t96 * t101 * 1.356336805555556E-7 + 1.0;
		F(17, 18) = t113 - b_wz * t63 * t78 * 2.0;
		F(17, 19) = t135 + t136;
		F(18, 11) = gy * (t106 + t119 - b_wx * t94 * 4.0)
				- gz
						* (t102 + t114 - b_wx * t63 * t80 * 2.0
								- b_wx * t67 * t78 * 2.0)
				- gx * (-t103 + t104 + t105 + b_wx * b_wy * t63 * t67 * 4.0);
		F(18, 12) = gy * (t110 + t126)
				+ gz * (t115 + t124 + t125 - b_wy * b_wz * t63 * t71 * 4.0)
				- gx * (-t97 + t107 + t108 + b_wx * b_wy * t63 * t71 * 4.0);
		F(18, 13) = gy * (t112 + t133 - b_wz * t94 * 4.0)
				+ gz * (t103 + t131 + t132 - b_wy * b_wz * t63 * t75 * 4.0)
				- gx
						* (-t102 + t111 + b_wz * t63 * t84 * 2.0
								+ b_wz * t75 * t78 * 2.0);
		F(18, 17) = t113 + b_wz * t63 * t78 * 2.0;
		F(18, 18) = t96 * t101 * (-1.356336805555556E-7)
				- t101 * t109 * 1.356336805555556E-7 + 1.0;
		F(18, 19) = t137 - b_wx * t63 * t78 * 2.0;
		F(19, 11) = gz * (t119 + t120 - b_wx * t94 * 4.0)
				+ gx * (t115 + t116 + t117 - t118)
				- gy
						* (-t102 + t114 + b_wx * t63 * t80 * 2.0
								+ b_wx * t67 * t78 * 2.0);
		F(19, 12) = gz * (t126 + t127 - b_wy * t94 * 4.0)
				- gx * (t102 - t121 - t122 + t123)
				- gy * (-t115 + t124 + t125 + b_wy * b_wz * t63 * t71 * 4.0);
		F(19, 13) = gz * (t133 + t134) + gx * (t97 + t128 + t129 - t130)
				- gy * (-t103 + t131 + t132 + b_wy * b_wz * t63 * t75 * 4.0);
		F(19, 17) = t135 - t136;
		F(19, 18) = t137 + b_wx * t63 * t78 * 2.0;
		F(19, 19) = t95 * t101 * (-1.356336805555556E-7)
				- t101 * t109 * 1.356336805555556E-7 + 1.0;
		F(20, 20) = 1.0;
		F(21, 21) = 1.0;
		F(22, 22) = 1.0;
		F(23, 23) = 1.0;
		F(24, 24) = 1.0;
		F(25, 25) = 1.0;

	}
}

