#pragma once
#ifndef CVL_H
#define CVL_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "source.h"
#include "area.h"

class COE
{
public:
	const int m = 4;
	const float mu = (4.0 * PI) * 1e-7f;
	const float epsilon = 8.85e-12f;
	const float eps_r = 1.f;
	const float eps0 = 8.85e-12f;
	const float PI = 3.141592653589793f;
	const float kappa_max = 8;
	const string filename = "coe.txt";

	//variables would not change after the set up
	float sigma_max;
	float num_layer;
	float d;
	float distance;
	float sigma;
	float alpha;
	float kappa;
	float cvl_coe;
	float c;

public:
	float set_distance(src, float xcoor);
	float set_sigma(src, float xcoor);
	float set_kappa(src, float xcoor);
	float set_alpha(src, float xcoor);
	float set_c(src, float xcoor);
	float set_coe(src, float xcoor);

public:
	COE(int, src);
	void checkout(src);
	~COE();
};

#endif