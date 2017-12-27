%computing the state transition jacobian for the loosely coupled visual
%odometry

syms x y z qw qx qy qz lambda b_dx b_dy b_dz b_wx b_wy b_wz b_ax b_ay b_az gx gy gz bias_accel_x bias_accel_y bias_accel_z bias_gyro_x bias_gyro_y bias_gyro_z dt

pos = [x;y;z]
quat = [qw;qx;qy;qz]
vel = [b_dx;b_dy;b_dy]
accel = [b_ax;b_ay;b_az]
omega = [b_wx;b_wy;b_wz]

pos_process = [pos + quaternionRotate(quat, lambda*(dt*vel + 0.5*dt*dt*accel))]

omega_norm = norm(omega);
v = omega/omega_norm;
theta = dt*omega_norm;
dq = [cos(theta/2); v*sin(theta/2)]
dq_inv = [cos(theta/2); -v*sin(theta/2)]

quat_process = multiplyQuaternions(quat, dq)