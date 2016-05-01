#pragma once
#ifndef CVL_H
#define CVL_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "source.h"
#include "E.h"
#include "H.h"

class E;
class H;

class cvl
{
public:
	//constants of equations
	const int m = 4;
	const float mu = (4.0 * PI) * 1e-7f;
	const float epsilon = 8.85e-12f;
	const float eps_r = 1.f;
	const float eps0 = 8.85e-12f;
	const float PI = 3.141592653589793f;
	const float kappa_max = 8;

	//variables would not change after the set up
	float sigma_max;
	int num_layer, side_sx, side_sy, ud_sx, ud_sy;
	float d;
	float *distance_full, *distance_half;
	float *sigma_full, *sigma_half;
	float *alpha_full, *alpha_half;
	float *kappa_full, *kappa_half;
	float *cvl_half_coe, *cvl_full_coe;

	//variables would change after every computation
	float *Hxzl, *Hxzr, *Hxzu, *Hxzd;
	float *Hyzl, *Hyzr, *Hyzu, *Hyzd;
	float *Exyl, *Exyr, *Exyu, *Exyd;
	float *Eyxl, *Eyxr, *Eyxu, *Eyxd;
	float *c_full, *c_half;
public:
	cvl(int , src);
	~cvl();

	void alloc(float* &, int, int);
	void set_distance(src);
	void set_sigma(src);
	void set_kappa(src);
	void set_alpha(src);
	void set_c(src);
	void set_half_coe(src);
	void set_full_coe(src);
	void checkout();
	void save2file();
	void save2file_e();
	void save2file_h();
	void save2file_coe();

	void cmp_cvlh(src, E, H, int);
	void cmp_cvle(src, E, H, int);
};

#endif