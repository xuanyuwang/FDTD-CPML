#pragma once
#ifndef CVL_H
#define CVL_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "source.h"

class cvl
{
public:
	const int m = 4;
	const float mu = (4.0 * PI) * 1e-7f;
	const float epsilon = 8.85e-12f;
	const float eps_r = 1.f;
	const float eps0 = 8.85e-12f;
	const float PI = 3.141592653589793f;
	const float kappa_max = 8;
	float sigma_max;
	int num_layer;
	float d;

	float *Hzxl, *Ezyl, *Exl, *Hyl;
	float *Hzxr, *Ezyr, *Exr, *Hyr;
	float *distance_E, *distance_H;
	float *sigma_E, *sigma_H;
	float *alpha_E, *alpha_H;
	float *kappa_E, *kappa_H;
	float *c_E, *c_H;
	float *cvl_H_coe, *cvl_E_coe;
public:
	cvl(int , src);
	~cvl();

	void set_HE(float*);
	void set_distance(src);
	void set_sigma(src);
	void set_kappa(src);
	void set_alpha(src);
	void set_c(src);
	void set_H_coe(src);
	void set_E_coe(src);
	void checkout();
	void save2file();
	void save2file_coe();

	void cmp_H_l(src source, float bd_E);
	void cmp_H_r(src source, float bd_E);
	void cmp_Hy(src s, float E_lbd, float E_rbd);

	void cmp_E_l(src);
	void cmp_E_r(src);
	void cmp_Ex(src s);
};

#endif