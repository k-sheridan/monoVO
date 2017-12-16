%computing the state transition jacobian for the loosely coupled visual
%odometry

syms x y z qw qx qy qz lambda b_dx b_dy b_dz b_wx b_wy b_wz b_ax b_ay b_az gx gy gz bias_accel_x bias_accel_y bias_accel_z bias_gyro_x bias_gyro_y bias_gyro_z dt

