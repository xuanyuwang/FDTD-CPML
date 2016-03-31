#pragma once
#ifndef SOURCE_H
#define SOURCE_H

#include <cmath>
#include <iostream>

using namespace std;

class src{
public:
	int size_space, size_time;
	float dt, dz;
	const float C = 3e8f;
public:
	src(int space, int time);
	void checkout();
	void cmp(int current_timestep, float* src_p);
};

#endif