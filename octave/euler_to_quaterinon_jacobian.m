syms r p y

t0 = cos(y * 0.5);
	 t1 = sin(y * 0.5);
	 t2 = cos(r * 0.5);
	 t3 = sin(r * 0.5);
	 t4 = cos(p * 0.5);
	 t5 = sin(p * 0.5);

	quat_fn = [t0 * t2 * t4 + t1 * t3 * t5;
	t0 * t3 * t4 - t1 * t2 * t5;
	t0 * t2 * t5 + t1 * t3 * t4;
	t1 * t2 * t4 - t0 * t3 * t5]

double(subs(quat_fn, [r, p, y], [pi/2,pi/2,0]))