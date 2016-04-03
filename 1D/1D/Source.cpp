#include "source.h"

src::src(int space, int time)
{
	dz = 0.015;
	dt = dz / (2 * C);
	size_space = space;
	size_time = time;
}

void src::checkout()
{
	cout << "dz = " << dz << endl;
	cout << "dt = " << dt << endl;
	cout << "size of space = " << size_space << endl;
	cout << "size of time = " << size_time << endl;
}

void src::cmp(int current_timestep, float* src_p)
{
	float T0, vt, val_src, time;

	time = current_timestep * dt;
	
	T0 = 3 * T;
	vt = (time - T0) / T;

	val_src = exp(-pow(vt, 2));

	*src_p = val_src;
}