#pragma once
#ifndef CVL_H
#define CVL_H

#include <iostream>
#include <cmath>
#include "source.h"
/*To Do*/
/*set_alpha*/
class cvl
{
public:
	const int m = 4;
	const float eps_r = 1.f;
	const float eps0 = 8.85e-12f;
	const float PI = 3.141592653589793f;
	const float sigma_max = (m + 1) / (sqrtf(eps_r) * 150 * PI);
	const float kappa_max = 8;
	float *Hzx, *Hxz, *Eyz, *Ezy;
	float dt, dz;
	int size_space;
	float sigma;
	float alpha;
	float d;
	float kappa;
	float c;
public:
	cvl(src source);
	~cvl();

	void set_sigma(float z, float z0);
	void set_kappa(float z, float z0);
	void set_alpha(float z, float z0);
	void set_c();
	//Hzx
	void cmp_H(int current_timestep, int z0, float* cvl_H, float *Ex);
	//Hxz
	void cmp_H(int current_timestep, int z0, float* cvl_H);
	//Eyz
	void cmp_E(int current_timestep, int z0, float *cvl_E, float *Hy);
	//Ezy
	void cmp_E(int current_timestep, int z0, float *cvl_E);
};

#endif