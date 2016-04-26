#pragma once
#ifndef SOURCE_H
#define SOURCE_H

#include <cmath>
#include <iostream>

using namespace std;

class src{
public:
	int num_grid, size_time;
	int size_x, size_y;
	float dt, dz;

	const float PI = 3.141592653589793f;
	const float C = 3e8f;
	const float T = 5e-10f;
	const float omega = (2 * PI) / T;
public:
	src(int space_x, int space_y, int time);
	void checkout();
	void cmp(int current_timestep, float* src_p);
};

#endif