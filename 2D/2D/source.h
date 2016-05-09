#pragma once
#ifndef SOURCE_H
#define SOURCE_H

#include <cmath>
#include <iostream>

using namespace std;

/************************************************************************/
/* Name: src
 * description:
 *	Member:
 *		num_grid: the number of Yee cells.
 *		size_time: the number of time steps
 *		size_x: the number of Yee cell of every row.
 *		size_y: the number of Yee cell of every colume.
 *		PI: constant.
 *		dt: the length of time step.
 *		dz: the size of ezch Yee cell.
 *		omega: the frequency of the source.
 *	Function: cmp(int, float*): assign the value of source of current
 *		time for a specific eletric field point.*/
/************************************************************************/
class src{
public:
	float num_grid, size_time;
	float size_x, size_y;
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