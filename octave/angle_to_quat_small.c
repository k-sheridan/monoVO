  t2 = dt*dt;
  t3 = fabs(b_wx);
  t4 = (b_wx/fabs(b_wx));
  t5 = fabs(b_wy);
  t6 = fabs(b_wz);
  t7 = t2*t2;
  t8 = (b_wy/fabs(b_wy));
  t9 = t3*t3;
  t10 = t5*t5;
  t11 = t6*t6;
  t12 = t9+t10+t11;
  t13 = (b_wz/fabs(b_wz));
  t14 = t2*t3*t4*(1.0/2.4E1);
  t22 = t3*t4*t7*t12*(1.0/9.6E2);
  t15 = t14-t22;
  t16 = t12*t12;
  t17 = t7*t16*2.604166666666667E-4;
  t18 = t2*t5*t8*(1.0/2.4E1);
  t23 = t5*t7*t8*t12*(1.0/9.6E2);
  t19 = t18-t23;
  t20 = t2*t6*t13*(1.0/2.4E1);
  t24 = t6*t7*t12*t13*(1.0/9.6E2);
  t21 = t20-t24;
  A0[0][0] = -t2*t3*t4+t3*t4*t7*t12*(1.0/9.6E1);
  A0[0][1] = -t2*t5*t8+t5*t7*t8*t12*(1.0/9.6E1);
  A0[0][2] = -t2*t6*t13+t6*t7*t12*t13*(1.0/9.6E1);
  A0[1][0] = t17-b_wx*t15-t2*t12*(1.0/4.8E1)+1.0/2.0;
  A0[1][1] = -b_wx*t19;
  A0[1][2] = -b_wx*t21;
  A0[2][0] = -b_wy*t15;
  A0[2][1] = t17-b_wy*t19-t2*t12*(1.0/4.8E1)+1.0/2.0;
  A0[2][2] = -b_wy*t21;
  A0[3][0] = -b_wz*t15;
  A0[3][1] = -b_wz*t19;
  A0[3][2] = t17-b_wz*t21-t2*t12*(1.0/4.8E1)+1.0/2.0;
