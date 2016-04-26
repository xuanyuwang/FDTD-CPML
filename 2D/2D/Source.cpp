#include "source.h"

src::src(int space_x, int space_y, int time)
{
	dz = 0.015;
	dt = dz / (2 * C);
	size_x = space_x;
	size_y = space_y;
	num_grid = size_x * size_y;
	size_time = time;
}

void src::checkout()
{
	cout << "dz = " << dz << endl;
	cout << "dt = " << dt << endl;
	cout << "space x-size = " << size_x << endl;
	cout << "space y-size = " << size_y << endl;
	cout << "grid number of space = " << num_grid << endl;
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