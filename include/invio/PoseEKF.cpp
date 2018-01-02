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
	this->state.setPosition(Eigen::Vector3d(0, 0, 0));
	this->state.setQuat(Eigen::Quaterniond(1, 0, 0, 0));
	this->state.setLambda(1.0);
	this->state.setBodyFrameVelocity(Eigen::Vector3d(0, 0, 0));
	this->state.setOmega(Eigen::Vector3d(0, 0, 0));
	this->state.setBodyFrameAcceleration(Eigen::Vector3d(0, 0, 0));
	this->state.setGravityVector(Eigen::Vector3d(0, 0, 9.8));

	//TODO set the initial biases as parameters
	this->state.setAccelerometerBias(Eigen::Vector3d(0, 0, 0));
	this->state.setGyroscopeBias(Eigen::Vector3d(0, 0, 0));

	//set the initial uncertainties
	this->state.Sigma(0, 0) = 10;
	this->state.Sigma(1, 1) = 10;
	this->state.Sigma(2, 2) = 10;

	this->state.Sigma(3, 3) = 10;
	this->state.Sigma(4, 4) = 10;
	this->state.Sigma(5, 5) = 10;
	this->state.Sigma(6, 6) = 10;

	this->state.Sigma(7, 7) = 0; //lambda is unknown when trying to estimate the scale

	this->state.Sigma(8, 8) = 10000;
	this->state.Sigma(9, 9) = 10000;
	this->state.Sigma(10, 10) = 10000;
	this->state.Sigma(11, 11) = 10000;
	this->state.Sigma(12, 12) = 10000;
	this->state.Sigma(13, 13) = 10000;
	this->state.Sigma(14, 14) = 10000;
	this->state.Sigma(15, 15) = 10000;
	this->state.Sigma(16, 16) = 10000;
	this->state.Sigma(17, 17) = 10000;
	this->state.Sigma(18, 18) = 10000;
	this->state.Sigma(19, 19) = 10000;

	this->state.Sigma(20, 20) = 10000;
	this->state.Sigma(21, 21) = 10000;
	this->state.Sigma(22, 22) = 10000;
	this->state.Sigma(23, 23) = 10000;
	this->state.Sigma(24, 24) = 10000;
	this->state.Sigma(25, 25) = 10000;

	ROS_DEBUG_STREAM("state: " << this->state.x);

}

Sophus::SE3d PoseEKF::getSE3(){
	Sophus::SE3d pose;

	pose = Sophus::SE3d(this->state.getQuat(), this->state.getPosition());

	return pose;
}

/*
 * this predicts the state forward in time and updates the uncertainties and correlations due to prediction uncertainty
 */
void PoseEKF::process(ros::Time to_time)
{
	double dt = (to_time - this->state.t).toSec();
	ROS_ASSERT(dt >= 0);

	ROS_DEBUG_STREAM("processing state: " << this->state.x << " with dt: " << dt);

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
	if (omega_norm > 1e-15)
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

	//use F to process the covariance matrix

	this->state.Sigma = F*this->state.Sigma*F.transpose() + this->computeProcessNoise(dt);

	ROS_DEBUG_STREAM("processed state: " << this->state.x << " with sigma: " << this->state.Sigma);

	this->state.t = to_time;
}

Eigen::Matrix<double, STATE_SIZE, STATE_SIZE> PoseEKF::computeProcessNoise(double dt)
{
	Eigen::Matrix<double, STATE_SIZE, STATE_SIZE> Q;
	Q.setZero(); // start with a fresh plate

	Q(0, 0) = 0.003*dt;
	Q(1, 1) = 0.003*dt;
	Q(2, 2) = 0.003*dt;

	Q(3, 3) = 0.003*dt;
	Q(4, 4) = 0.003*dt;
	Q(5, 5) = 0.003*dt;
	Q(6, 6) = 0.003*dt;

	Q(7, 7) = 0*dt; //

	Q(8, 8) = 2*dt;
	Q(9, 9) = 2*dt;
	Q(10, 10) = 2*dt;
	Q(11, 11) = 2*dt;
	Q(12, 12) = 2*dt;
	Q(13, 13) = 2*dt;

	Q(14, 14) = 4*dt;
	Q(15, 15) = 4*dt;
	Q(16, 16) = 4*dt;

	Q(17, 17) = 0.003*dt;
	Q(18, 18) = 0.003*dt;
	Q(19, 19) = 0.003*dt;

	Q(20, 20) = 0.005*dt;
	Q(21, 21) = 0.005*dt;
	Q(22, 22) = 0.005*dt;
	Q(23, 23) = 0.005*dt;
	Q(24, 24) = 0.005*dt;
	Q(25, 25) = 0.005*dt;

	return Q;
}

/*
 * need to convert the euler covariance matrix to a quaternion covariance matrix before update can
 * be performed. This is a completely linear process other wise
 */
void PoseEKF::updateWithVOPose(Sophus::SE3d pose,
		Eigen::Matrix<double, 6, 6> cov, ros::Time t_measured)
{
	ROS_DEBUG("updating pose");

	//we measure an euler angle, but our state stores its orientation in a quaternion

	Eigen::Matrix<double, 7, 6> J;
	J.setZero();
	this->computeAngleAxis2QuaternionSpaceJacobian(pose.log(), J);

	Eigen::Matrix<double, 7, 7> transformed_cov = J*cov*J.transpose(); // finally get the covariance in terms of quaternions

	ROS_DEBUG_STREAM("updating with measurement: " << pose.log() << " and cov: " << cov);
	ROS_DEBUG_STREAM("transformed cov: " << transformed_cov);


	Eigen::Matrix<double, 7, 1> residual;
	residual << pose.translation().x() - this->state.x(0), pose.translation().y() - this->state.x(1), pose.translation().z() - this->state.x(2), pose.unit_quaternion().w() - this->state.x(3),
			pose.unit_quaternion().x() - this->state.x(4), pose.unit_quaternion().y() - this->state.x(5), pose.unit_quaternion().z() - this->state.x(6);

	Eigen::Matrix<double, 7, STATE_SIZE> H;
	H.setIdentity(); // the measurement transition is identity

	/*
	 * y = measurement - expected_measurement;
    S = R + H*Sigma*H';
    K = Sigma*H'/(S);

    mu = mu + K*y;

    % use Josephs form: P = ( I ? KH) P (I ? KH)' + KRK'
    I_KH = (eye(size(Sigma)) - K*H);
    Sigma = I_KH*Sigma*I_KH' + K*R*K'
	 */

	Eigen::Matrix<double, 7, 7> S = transformed_cov + H*this->state.Sigma*H.transpose();
	//x*A = b is equivalent to A.transpose() * z = b.transpose(); x = z.transpose()


	Eigen::Matrix<double, STATE_SIZE, 7> K = this->state.Sigma * H.transpose() * S.lu().inverse(); // TODO check efficiency

	this->state.x = this->state.x + K*residual;

	//normalize the quaternion again due to the linearization.
	this->state.setQuat(this->state.getQuat().normalized());


	// update the covariance
	Eigen::Matrix<double, STATE_SIZE, STATE_SIZE> I_KH = Eigen::MatrixXd::Identity(STATE_SIZE, STATE_SIZE) - K * H;

	this->state.Sigma = I_KH*this->state.Sigma*I_KH.transpose() + K*transformed_cov*K.transpose();

}

void PoseEKF::computeAngleAxis2QuaternionSpaceJacobian(Eigen::Matrix<double, 6, 1> x,
		Eigen::Matrix<double, 7, 6>& J)
{
	// this is hand computed
	//x = x
	//y = y
	//z = z
	//qw = cos(sqrt(r^2+p^2+y^2)/2)
	//qx = (r/sqrt(r^2+p^2+y^2)) * sin(sqrt(r^2+p^2+y^2)/2)
	//qy = (p/sqrt(r^2+p^2+y^2)) * sin(sqrt(r^2+p^2+y^2)/2)
	//qz = (y/sqrt(r^2+p^2+y^2)) * sin(sqrt(r^2+p^2+y^2)/2)

	//jacobian from here: https://math.stackexchange.com/questions/1996586/transform-se3-pose-covariance

	Eigen::Vector3d r = Eigen::Vector3d(x(3),x(4),x(5));

	J(0, 0) = 1;
	J(1, 1) = 1;
	J(2, 2) = 1;

	double phi = r.norm();

	if(phi > 0.00000000000001)
	{
		double alpha = sin(phi/2)/phi;

		J.block<3, 3>(3, 3) = ((0.5 * cos(phi/2) - alpha)/(phi*phi)) * r * r.transpose() + alpha*Eigen::MatrixXd::Identity(3, 3);

		J.block<1, 3>(6, 3) = -alpha/2 * r.transpose();
	}
	else
	{
		double alpha = 0.5;

		J.block<3, 3>(3, 3) = -0.04167 * r * r.transpose() + alpha*Eigen::MatrixXd::Identity(3, 3);

		J.block<1, 3>(6, 3) = -alpha/2 * r.transpose();
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
		double t2 = dt*dt;
		double t3 = b_dy*dt;
		double   t4 = b_ay*t2*(1.0/2.0);
		double   t5 = t3+t4;
		double   t6 = b_dz*dt;
		double   t7 = b_az*t2*(1.0/2.0);
		double   t8 = t6+t7;
		double   t9 = b_dx*dt;
		double   t10 = b_ax*t2*(1.0/2.0);
		double   t11 = t9+t10;
		double   t12 = qy*qy;
		double   t13 = t12*2.0;
		double   t14 = qz*qz;
		double   t15 = t14*2.0;
		double   t16 = t13+t15-1.0;
		double 	  t17 = qw*qz*2.0;
		double 	  t18 = qx*qy*2.0;
		double 	  t19 = t17+t18;
		double 	  t20 = qw*qy;
		double 	  t54 = qx*qz;
		double 	  t21 = t20-t54;
		double 	  t22 = lambda*qx*t8*2.0;
		double 	  t23 = lambda*qz*t8*2.0;
		double 	  t24 = b_dy*2.0;
		double 	  t25 = b_ay*dt;
		double 	  t26 = t24+t25;
		double 	  t27 = b_dx*2.0;
		double 	  t28 = b_ax*dt;
		double 	  t29 = t27+t28;
		double 	  t30 = b_dz*2.0;
		double 	  t31 = b_az*dt;
		double 	  t32 = t30+t31;
		double 	  t33 = qw*qz;
		double 	  t34 = qx*qy;
		double 	  t35 = qx*qx;
		double 	  t36 = t35*2.0;
		double 	  t37 = t15+t36-1.0;
		double 	  t38 = qw*qx*2.0;
		double 	  t39 = qy*qz*2.0;
		double 	  t40 = t38+t39;
		double 	  t41 = t33-t34;
		double 	  t42 = lambda*qx*t5*2.0;
		double 	  t43 = lambda*qy*t11*2.0;
		double 	  t44 = lambda*qw*t5*2.0;
		double 	  t45 = lambda*qz*t5*2.0;
		double 	  t46 = lambda*qx*t11*2.0;
		double 	  t47 = lambda*qy*t5*2.0;
		double 	  t48 = qw*qy*2.0;
		double 	  t49 = qx*qz*2.0;
		double 	  t50 = t48+t49;
		double 	  t51 = qw*qx;
		double 	  t52 = qy*qz;
		double 	  t53 = t13+t36-1.0;
		double 		  t55 = t51-t52;
		double 	  t56 = fabs(b_wx);
		double 	  t57 = fabs(b_wy);
		double 	  t58 = fabs(b_wz);
		double 	  t59 = t56*t56;
		double 	  t60 = t57*t57;
		double 	  t61 = t58*t58;
		double 	  t62 = t59+t60+t61;
		double 	  t63 = sqrt(t62);
		double 	  t64 = dt*t63*(1.0/2.0);
		double 	  t65 = sin(t64);
		double 	  t66 = 1.0/sqrt(t62);
		double 	  t67 = (b_wx/fabs(b_wx));
		double   t68 = 1.0/pow(t62,3.0/2.0);
		double 	  t69 = cos(t64);
		double 	  t70 = 1.0/t62;
		double 	  t71 = (b_wy/fabs(b_wy));
		double 	  t72 = (b_wz/fabs(b_wz));
		double 	  t73 = b_wz*t65*t66;
		double 	  t74 = b_wx*t65*t66;
		double 	  t75 = qw*t65*t66;
		double 	  t76 = b_wy*t65*t66;
		double 	  t77 = qy*t65*t66;
		double 	  t78 = t65*t65;
		double 	  t79 = 1.0/(t62*t62);
		double 	  t80 = b_wy*b_wy;
		double 	  t81 = b_wz*b_wz;
		double 	  t82 = t69*t69;
		double 	  t83 = dt*t63;
		double 	  t84 = sin(t83);
		double 	  t85 = t66*t84;
		double 	  t86 = b_wx*t70*t78*2.0;
		double 	  t87 = t56*t67*t78*t79*t81*4.0;
		double 	  t88 = b_wx*b_wx;
		double 	  t89 = b_wy*t70*t78*2.0;
		double 	  t90 = b_wz*t56*t65*t67*t68*t69*2.0;
		double 		  t91 = b_wz*dt*t56*t67*t70*t78;
		double 		  t92 = b_wx*b_wy*dt*t56*t65*t67*t68*t69*2.0;
		double 		  t93 = dt*t57*t65*t68*t69*t71*t81*2.0;
		double 		  t94 = b_wz*t57*t65*t68*t69*t71*2.0;
		double 		  t95 = b_wz*dt*t57*t70*t71*t78;
		double 		  t96 = b_wx*b_wy*dt*t57*t65*t68*t69*t71*2.0;
		double 		  t97 = b_wz*t70*t78*2.0;
		double 	  t98 = b_wz*t70*t78*4.0;
		double 	  t99 = dt*t58*t65*t68*t69*t72*t81*2.0;
		double 	  t100 = b_wz*dt*t58*t70*t72*t82;
		double 	  t101 = b_wx*b_wy*t58*t72*t78*t79*4.0;
		double 	  t102 = b_wx*b_wy*t70*t78*2.0;
		double 	  t103 = b_wx*t70*t78*4.0;
		double 	  t104 = t56*t67*t78*t79*t80*4.0;
		double 	  t105 = dt*t56*t65*t67*t68*t69*t88*2.0;
		double 	  t106 = b_wx*dt*t56*t67*t70*t82;
		double 	  t107 = b_wy*b_wz*t56*t67*t78*t79*4.0;
		double 	  t108 = b_wy*dt*t56*t67*t70*t82;
		double 	  t109 = b_wx*b_wz*dt*t56*t65*t67*t68*t69*2.0;
		double 	  t110 = b_wy*t70*t78*4.0;
		double 	  t111 = dt*t57*t65*t68*t69*t71*t88*2.0;
		double 	  t112 = dt*t57*t65*t68*t69*t71*t80*2.0;
		double 	  t113 = b_wy*dt*t57*t70*t71*t82;
		double 	  t114 = b_wx*b_wz*dt*t57*t65*t68*t69*t71*2.0;
		double 	  t115 = b_wx*t57*t65*t68*t69*t71*2.0;
		double 	  t116 = b_wx*dt*t57*t70*t71*t78;
		double 	  t117 = b_wy*b_wz*dt*t57*t65*t68*t69*t71*2.0;
		double 	  t118 = dt*t58*t65*t68*t69*t72*t88*2.0;
		double 	  t119 = dt*t58*t65*t68*t69*t72*t80*2.0;
		double 	  t120 = b_wy*dt*t58*t70*t72*t82;
		double 	  t121 = b_wx*b_wz*dt*t58*t65*t68*t69*t72*2.0;
		double 	  t122 = b_wx*t58*t65*t68*t69*t72*2.0;
		double 	  t123 = b_wx*dt*t58*t70*t72*t78;
		double 	  t124 = b_wy*b_wz*dt*t58*t65*t68*t69*t72*2.0;
		double 	  t125 = b_wy*t66*t84;
		double 	  t126 = b_wx*b_wz*t70*t78*2.0;
		double   t127 = b_wy*b_wz*t70*t78*2.0;
		F(0,0) = 1.0;
		F(0,3) = t45-lambda*qy*t8*2.0;
		F(0,4) = t23+t47;
		F(0,5) = t42-lambda*qw*t8*2.0-lambda*qy*t11*4.0;
		F(0,6) = t22+t44-lambda*qz*t11*4.0;
		F(0,7) = dt*t32*(t48-qx*qz*2.0)*(-1.0/2.0)-dt*t16*t29*(1.0/2.0)+dt*t19*t26*(1.0/2.0);
		F(0,8) = -dt*lambda*t16;
		F(0,9) = dt*lambda*t19;
		F(0,10) = dt*lambda*t21*-2.0;
		F(0,14) = lambda*t2*t16*(-1.0/2.0);
		F(0,15) = lambda*t2*(t33+t34);
		F(0,16) = -lambda*t2*t21;
		F(1,1) = 1.0;
		F(1,3) = t22-lambda*qz*t11*2.0;
		F(1,4) = t43+lambda*qw*t8*2.0-lambda*qx*t5*4.0;
		F(1,5) = t23+t46;
		F(1,6) = lambda*qw*t11*-2.0+lambda*qy*t8*2.0-lambda*qz*t5*4.0;
		F(1,7) = dt*t26*t37*(-1.0/2.0)+dt*t32*t40*(1.0/2.0)-dt*t29*(t17-t18)*(1.0/2.0);
		F(1,8) = dt*lambda*t41*-2.0;
		F(1,9) = -dt*lambda*t37;
		F(1,10) = dt*lambda*t40;
		F(1,14) = -lambda*t2*t41;
		F(1,15) = lambda*t2*t37*(-1.0/2.0);
		F(1,16) = lambda*t2*(t51+t52);
		F(2,2) = 1.0;
		F(2,3) = -t42+t43;
		F(2,4) = -t44-lambda*qx*t8*4.0+lambda*qz*t11*2.0;
		F(2,5) = t45+lambda*qw*t11*2.0-lambda*qy*t8*4.0;
		F(2,6) = t46+t47;
		F(2,7) = dt*t29*t50*(1.0/2.0)-dt*t32*t53*(1.0/2.0)-dt*t26*(t38-t39)*(1.0/2.0);
		F(2,8) = dt*lambda*t50;
		F(2,9) = dt*lambda*t55*-2.0;
		F(2,10) = -dt*lambda*t53;
		F(2,14) = lambda*t2*(t20+t54);
		F(2,15) = -lambda*t2*t55;
		F(2,16) = lambda*t2*t53*(-1.0/2.0);
		F(3,3) = t69;
		F(3,4) = -b_wx*t65*t66;
		F(3,5) = -b_wy*t65*t66;
		F(3,6) = -b_wz*t65*t66;
		F(3,11) = -qx*t65*t66+b_wx*qx*t56*t65*t67*t68+b_wy*qy*t56*t65*t67*t68+b_wz*qz*t56*t65*t67*t68-dt*qw*t56*t65*t66*t67*(1.0/2.0)-b_wx*dt*qx*t56*t67*t69*t70*(1.0/2.0)-b_wy*dt*qy*t56*t67*t69*t70*(1.0/2.0)-b_wz*dt*qz*t56*t67*t69*t70*(1.0/2.0);
		F(3,12) = -qy*t65*t66+b_wx*qx*t57*t65*t68*t71+b_wy*qy*t57*t65*t68*t71+b_wz*qz*t57*t65*t68*t71-dt*qw*t57*t65*t66*t71*(1.0/2.0)-b_wx*dt*qx*t57*t69*t70*t71*(1.0/2.0)-b_wy*dt*qy*t57*t69*t70*t71*(1.0/2.0)-b_wz*dt*qz*t57*t69*t70*t71*(1.0/2.0);
		F(3,13) = -qz*t65*t66+b_wx*qx*t58*t65*t68*t72+b_wy*qy*t58*t65*t68*t72+b_wz*qz*t58*t65*t68*t72-dt*qw*t58*t65*t66*t72*(1.0/2.0)-b_wx*dt*qx*t58*t69*t70*t72*(1.0/2.0)-b_wy*dt*qy*t58*t69*t70*t72*(1.0/2.0)-b_wz*dt*qz*t58*t69*t70*t72*(1.0/2.0);
		F(4,3) = t74;
		F(4,4) = t69;
		F(4,5) = t73;
		F(4,6) = -b_wy*t65*t66;
		F(4,11) = t75-b_wx*qw*t56*t65*t67*t68-b_wz*qy*t56*t65*t67*t68+b_wy*qz*t56*t65*t67*t68-dt*qx*t56*t65*t66*t67*(1.0/2.0)+b_wx*dt*qw*t56*t67*t69*t70*(1.0/2.0)+b_wz*dt*qy*t56*t67*t69*t70*(1.0/2.0)-b_wy*dt*qz*t56*t67*t69*t70*(1.0/2.0);
		F(4,12) = -qz*t65*t66-b_wx*qw*t57*t65*t68*t71-b_wz*qy*t57*t65*t68*t71+b_wy*qz*t57*t65*t68*t71-dt*qx*t57*t65*t66*t71*(1.0/2.0)+b_wx*dt*qw*t57*t69*t70*t71*(1.0/2.0)+b_wz*dt*qy*t57*t69*t70*t71*(1.0/2.0)-b_wy*dt*qz*t57*t69*t70*t71*(1.0/2.0);
		F(4,13) = t77-b_wx*qw*t58*t65*t68*t72-b_wz*qy*t58*t65*t68*t72+b_wy*qz*t58*t65*t68*t72-dt*qx*t58*t65*t66*t72*(1.0/2.0)+b_wx*dt*qw*t58*t69*t70*t72*(1.0/2.0)+b_wz*dt*qy*t58*t69*t70*t72*(1.0/2.0)-b_wy*dt*qz*t58*t69*t70*t72*(1.0/2.0);
		F(5,3) = t76;
		F(5,4) = -t73;
		F(5,5) = t69;
		F(5,6) = t74;
		F(5,11) = qz*t65*t66-b_wy*qw*t56*t65*t67*t68+b_wz*qx*t56*t65*t67*t68-b_wx*qz*t56*t65*t67*t68-dt*qy*t56*t65*t66*t67*(1.0/2.0)+b_wy*dt*qw*t56*t67*t69*t70*(1.0/2.0)-b_wz*dt*qx*t56*t67*t69*t70*(1.0/2.0)+b_wx*dt*qz*t56*t67*t69*t70*(1.0/2.0);
		F(5,12) = t75-b_wy*qw*t57*t65*t68*t71+b_wz*qx*t57*t65*t68*t71-b_wx*qz*t57*t65*t68*t71-dt*qy*t57*t65*t66*t71*(1.0/2.0)+b_wy*dt*qw*t57*t69*t70*t71*(1.0/2.0)-b_wz*dt*qx*t57*t69*t70*t71*(1.0/2.0)+b_wx*dt*qz*t57*t69*t70*t71*(1.0/2.0);
		F(5,13) = -qx*t65*t66-b_wy*qw*t58*t65*t68*t72+b_wz*qx*t58*t65*t68*t72-b_wx*qz*t58*t65*t68*t72-dt*qy*t58*t65*t66*t72*(1.0/2.0)+b_wy*dt*qw*t58*t69*t70*t72*(1.0/2.0)-b_wz*dt*qx*t58*t69*t70*t72*(1.0/2.0)+b_wx*dt*qz*t58*t69*t70*t72*(1.0/2.0);
		F(6,3) = t73;
		F(6,4) = t76;
		F(6,5) = -t74;
		F(6,6) = t69;
		F(6,11) = -t77-b_wz*qw*t56*t65*t67*t68-b_wy*qx*t56*t65*t67*t68+b_wx*qy*t56*t65*t67*t68-dt*qz*t56*t65*t66*t67*(1.0/2.0)+b_wz*dt*qw*t56*t67*t69*t70*(1.0/2.0)+b_wy*dt*qx*t56*t67*t69*t70*(1.0/2.0)-b_wx*dt*qy*t56*t67*t69*t70*(1.0/2.0);
		F(6,12) = qx*t65*t66-b_wz*qw*t57*t65*t68*t71-b_wy*qx*t57*t65*t68*t71+b_wx*qy*t57*t65*t68*t71-dt*qz*t57*t65*t66*t71*(1.0/2.0)+b_wz*dt*qw*t57*t69*t70*t71*(1.0/2.0)+b_wy*dt*qx*t57*t69*t70*t71*(1.0/2.0)-b_wx*dt*qy*t57*t69*t70*t71*(1.0/2.0);
		F(6,13) = t75-b_wz*qw*t58*t65*t68*t72-b_wy*qx*t58*t65*t68*t72+b_wx*qy*t58*t65*t68*t72-dt*qz*t58*t65*t66*t72*(1.0/2.0)+b_wz*dt*qw*t58*t69*t70*t72*(1.0/2.0)+b_wy*dt*qx*t58*t69*t70*t72*(1.0/2.0)-b_wx*dt*qy*t58*t69*t70*t72*(1.0/2.0);
		F(7,7) = 1.0;
		F(8,8) = 1.0;
		F(8,14) = dt;
		F(9,9) = 1.0;
		F(9,15) = dt;
		F(10,10) = 1.0;
		F(10,16) = dt;
		F(11,11) = 1.0;
		F(12,12) = 1.0;
		F(13,13) = 1.0;
		F(14,14) = 1.0;
		F(15,15) = 1.0;
		F(16,16) = 1.0;
		F(17,11) = gx*(t87+t104-dt*t56*t65*t67*t68*t69*t80*2.0-dt*t56*t65*t67*t68*t69*t81*2.0)+gz*(t97+t108+t109-b_wy*t56*t65*t67*t68*t69*2.0-b_wx*b_wz*t56*t67*t78*t79*4.0-b_wy*dt*t56*t67*t70*t78)+gy*(t89+t90+t91+t92-b_wx*b_wy*t56*t67*t78*t79*4.0-b_wz*dt*t56*t67*t70*t82);
		F(17,12) = gz*(t85+t113+t114-b_wy*t57*t65*t68*t69*t71*2.0-b_wx*b_wz*t57*t71*t78*t79*4.0-b_wy*dt*t57*t70*t71*t78)-gx*(t93+t110+t112-t57*t71*t78*t79*t80*4.0-t57*t71*t78*t79*t81*4.0)+gy*(t86+t94+t95+t96-b_wx*b_wy*t57*t71*t78*t79*4.0-b_wz*dt*t57*t70*t71*t82);
		F(17,13) = -gy*(t85+t100+t101-b_wz*t58*t65*t68*t69*t72*2.0-b_wz*dt*t58*t70*t72*t78-b_wx*b_wy*dt*t58*t65*t68*t69*t72*2.0)+gz*(t86+t120+t121-b_wy*t58*t65*t68*t69*t72*2.0-b_wx*b_wz*t58*t72*t78*t79*4.0-b_wy*dt*t58*t70*t72*t78)-gx*(t98+t99+t119-t58*t72*t78*t79*t80*4.0-t58*t72*t78*t79*t81*4.0);
		F(17,17) = t70*t78*t80*-2.0-t70*t78*t81*2.0+1.0;
		F(17,18) = t102-b_wz*t66*t84;
		F(17,19) = t125+t126;
		F(18,11) = -gy*(-t87+t103+t105-t56*t67*t78*t79*t88*4.0+dt*t56*t65*t67*t68*t69*t81*2.0)-gz*(t85+t106+t107-b_wx*t56*t65*t67*t68*t69*2.0-b_wx*dt*t56*t67*t70*t78-b_wy*b_wz*dt*t56*t65*t67*t68*t69*2.0)+gx*(t89-t90-t91+t92-b_wx*b_wy*t56*t67*t78*t79*4.0+b_wz*dt*t56*t67*t70*t82);
		F(18,12) = -gy*(t93+t111-t57*t71*t78*t79*t81*4.0-t57*t71*t78*t79*t88*4.0)+gx*(t86-t94-t95+t96-b_wx*b_wy*t57*t71*t78*t79*4.0+b_wz*dt*t57*t70*t71*t82)+gz*(t97+t115+t116+t117-b_wy*b_wz*t57*t71*t78*t79*4.0-b_wx*dt*t57*t70*t71*t82);
		F(18,13) = gx*(t85+t100-t101-b_wz*t58*t65*t68*t69*t72*2.0-b_wz*dt*t58*t70*t72*t78+b_wx*b_wy*dt*t58*t65*t68*t69*t72*2.0)-gy*(t98+t99+t118-t58*t72*t78*t79*t81*4.0-t58*t72*t78*t79*t88*4.0)+gz*(t89+t122+t123+t124-b_wy*b_wz*t58*t72*t78*t79*4.0-b_wx*dt*t58*t70*t72*t82);
		F(18,17) = t102+b_wz*t66*t84;
		F(18,18) = t70*t78*t81*-2.0-t70*t78*t88*2.0+1.0;
		F(18,19) = t127-b_wx*t66*t84;
		F(19,11) = -gz*(t103-t104+t105-t56*t67*t78*t79*t88*4.0+dt*t56*t65*t67*t68*t69*t80*2.0)+gy*(t85+t106-t107-b_wx*t56*t65*t67*t68*t69*2.0-b_wx*dt*t56*t67*t70*t78+b_wy*b_wz*dt*t56*t65*t67*t68*t69*2.0)+gx*(t97-t108+t109+b_wy*t56*t65*t67*t68*t69*2.0-b_wx*b_wz*t56*t67*t78*t79*4.0+b_wy*dt*t56*t67*t70*t78);
		F(19,12) = gy*(t97-t115-t116+t117-b_wy*b_wz*t57*t71*t78*t79*4.0+b_wx*dt*t57*t70*t71*t82)-gz*(t110+t111+t112-t57*t71*t78*t79*t80*4.0-t57*t71*t78*t79*t88*4.0)-gx*(t85+t113-t114-b_wy*t57*t65*t68*t69*t71*2.0+b_wx*b_wz*t57*t71*t78*t79*4.0-b_wy*dt*t57*t70*t71*t78);
		F(19,13) = -gz*(t118+t119-t58*t72*t78*t79*t80*4.0-t58*t72*t78*t79*t88*4.0)+gy*(t89-t122-t123+t124-b_wy*b_wz*t58*t72*t78*t79*4.0+b_wx*dt*t58*t70*t72*t82)+gx*(t86-t120+t121+b_wy*t58*t65*t68*t69*t72*2.0-b_wx*b_wz*t58*t72*t78*t79*4.0+b_wy*dt*t58*t70*t72*t78);
		F(19,17) = -t125+t126;
		F(19,18) = t127+b_wx*t66*t84;
		F(19,19) = t70*t78*t80*-2.0-t70*t78*t88*2.0+1.0;
		F(20,20) = 1.0;
		F(21,21) = 1.0;
		F(22,22) = 1.0;
		F(23,23) = 1.0;
		F(24,24) = 1.0;
		F(25,25) = 1.0;


	}
	else
	{
		double t2 = dt*dt;
		double t3 = b_dy*dt;
		double t4 = b_ay*t2*(1.0/2.0);
		double t5 = t3+t4;
		double t6 = b_dz*dt;
		double t7 = b_az*t2*(1.0/2.0);
		double t8 = t6+t7;
		double t9 = b_dx*dt;
		double t10 = b_ax*t2*(1.0/2.0);
		double t11 = t9+t10;
		double t12 = qy*qy;
		double t13 = t12*2.0;
		double t14 = qz*qz;
		double 	  t15 = t14*2.0;
		double 	  t16 = t13+t15-1.0;
		double 	  t17 = qw*qz*2.0;
		double   t18 = qx*qy*2.0;
		double   t19 = t17+t18;
		double   t20 = qw*qy;
		double   t54 = qx*qz;
		double   t21 = t20-t54;
		double   t22 = lambda*qx*t8*2.0;
		double   t23 = lambda*qz*t8*2.0;
		double   t24 = b_dy*2.0;
		double  t25 = b_ay*dt;
		double  t26 = t24+t25;
		double  t27 = b_dx*2.0;
		double  t28 = b_ax*dt;
		double  t29 = t27+t28;
		double   t30 = b_dz*2.0;
		double  t31 = b_az*dt;
		double  t32 = t30+t31;
		double  t33 = qw*qz;
		double   t34 = qx*qy;
		double   t35 = qx*qx;
		double  t36 = t35*2.0;
		double   t37 = t15+t36-1.0;
		double  t38 = qw*qx*2.0;
		double  t39 = qy*qz*2.0;
		double  t40 = t38+t39;
		double  t41 = t33-t34;
		double  t42 = lambda*qx*t5*2.0;
		double  t43 = lambda*qy*t11*2.0;
		double  t44 = lambda*qw*t5*2.0;
		double   t45 = lambda*qz*t5*2.0;
		double   t46 = lambda*qx*t11*2.0;
		double  t47 = lambda*qy*t5*2.0;
		double  t48 = qw*qy*2.0;
		double   t49 = qx*qz*2.0;
		double  t50 = t48+t49;
		double  t51 = qw*qx;
		double  t52 = qy*qz;
		double  t53 = t13+t36-1.0;
		double  t55 = t51-t52;
		F(0,0) = 1.0;
		F(0,3) = t45-lambda*qy*t8*2.0;
		F(0,4) = t23+t47;
		F(0,5) = t42-lambda*qw*t8*2.0-lambda*qy*t11*4.0;
		F(0,6) = t22+t44-lambda*qz*t11*4.0;
		F(0,7) = dt*t32*(t48-qx*qz*2.0)*(-1.0/2.0)-dt*t16*t29*(1.0/2.0)+dt*t19*t26*(1.0/2.0);
		F(0,8) = -dt*lambda*t16;
		F(0,9) = dt*lambda*t19;
		F(0,10) = dt*lambda*t21*-2.0;
		F(0,14) = lambda*t2*t16*(-1.0/2.0);
		F(0,15) = lambda*t2*(t33+t34);
		F(0,16) = -lambda*t2*t21;
		F(1,1) = 1.0;
		F(1,3) = t22-lambda*qz*t11*2.0;
		F(1,4) = t43+lambda*qw*t8*2.0-lambda*qx*t5*4.0;
		F(1,5) = t23+t46;
		F(1,6) = lambda*qw*t11*-2.0+lambda*qy*t8*2.0-lambda*qz*t5*4.0;
		F(1,7) = dt*t26*t37*(-1.0/2.0)+dt*t32*t40*(1.0/2.0)-dt*t29*(t17-t18)*(1.0/2.0);
		F(1,8) = dt*lambda*t41*-2.0;
		F(1,9) = -dt*lambda*t37;
		F(1,10) = dt*lambda*t40;
		F(1,14) = -lambda*t2*t41;
		F(1,15) = lambda*t2*t37*(-1.0/2.0);
		F(1,16) = lambda*t2*(t51+t52);
		F(2,2) = 1.0;
		F(2,3) = -t42+t43;
		F(2,4) = -t44-lambda*qx*t8*4.0+lambda*qz*t11*2.0;
		F(2,5) = t45+lambda*qw*t11*2.0-lambda*qy*t8*4.0;
		F(2,6) = t46+t47;
		F(2,7) = dt*t29*t50*(1.0/2.0)-dt*t32*t53*(1.0/2.0)-dt*t26*(t38-t39)*(1.0/2.0);
		F(2,8) = dt*lambda*t50;
		F(2,9) = dt*lambda*t55*-2.0;
		F(2,10) = -dt*lambda*t53;
		F(2,14) = lambda*t2*(t20+t54);
		F(2,15) = -lambda*t2*t55;
		F(2,16) = lambda*t2*t53*(-1.0/2.0);
		F(3,3) = 1.0;
		F(4,4) = 1.0;
		F(5,5) = 1.0;
		F(6,6) = 1.0;
		F(7,7) = 1.0;
		F(8,8) = 1.0;
		F(8,14) = dt;
		F(9,9) = 1.0;
		F(9,15) = dt;
		F(10,10) = 1.0;
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

