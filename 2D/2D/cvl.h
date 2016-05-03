#pragma once
#ifndef CVL_H
#define CVL_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "source.h"
#include "E.h"
#include "H.h"
#include "TOOLS.h"

class E;
class HX;
class HY;

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

	//variables would not change after the set up
	float sigma_max;
	int num_layer, side_sx, side_sy, ud_sx, ud_sy;
	float d;
	float *distance_full, *distance_half;
	float *sigma_full, *sigma_half;
	float *alpha_full, *alpha_half;
	float *kappa_full, *kappa_half;
	float *cvl_half_coe, *cvl_full_coe;
	float *c_full, *c_half;

public:
	void set_distance(src);
	void set_sigma(src);
	void set_kappa(src);
	void set_alpha(src);
	void set_c(src);
	void set_half_coe(src);
	void set_full_coe(src);

public:
	COE(int, src);
	void checkout();
	~COE();
};

class HXZL
{
public:
	float *p;
	int width, height;
public:
	HXZL(COE, src s);
	void cmp(COE, src, E, int);
	void save2file();
};

class HXZR
{
public:
	float *p;
	int width, height;
public:
	HXZR(COE, src s);
	void cmp(COE, src, E, int);
	void save2file();
};

class HXZU
{
public:
	float *p;
	int width, height;
public:
	HXZU(COE, src s);
	void cmp(COE, src, E, int);
	void save2file();
};

class HXZD
{
public:
	float *p;
	int width, height;
public:
	HXZD(COE, src s);
	void cmp(COE, src, E, int);
	void save2file();
};

class HYZL
{
public:
	float *p;
	int width, height;
public:
	HYZL(COE, src s);
	void cmp(COE, src, E, int);
	void save2file();
};

class HYZR
{
public:
	float *p;
	int width, height;
public:
	HYZR(COE, src s);
	void cmp(COE, src, E, int);
	void save2file();
};

class HYZU
{
public:
	float *p;
	int width, height;
public:
	HYZU(COE, src s);
	void cmp(COE, src, E, int);
	void save2file();
};

class HYZD
{
public:
	float *p;
	int width, height;
public:
	HYZD(COE, src s);
	void cmp(COE, src, E, int);
	void save2file();
};

class EXYL
{
public:
	float *p;
	int width, height;
public:
	EXYL(COE, src s);
	void cmp(COE, src, HY, int);
	void save2file();
};

class EXYR
{
public:
	float *p;
	int width, height;
public:
	EXYR(COE, src s);
	void cmp(COE, src, HY, int);
	void save2file();
};

class EXYU
{
public:
	float *p;
	int width, height;
public:
	EXYU(COE, src s);
	void cmp(COE, src, HY, int);
	void save2file();
};

class EXYD
{
public:
	float *p;
	int width, height;
public:
	EXYD(COE, src s);
	void cmp(COE, src, HY, int);
	void save2file();
};

class EYXL
{
public:
	float *p;
	int width, height;
public:
	EYXL(COE, src s);
	void cmp(COE, src, HX, int);
	void save2file();
};

class EYXR
{
public:
	float *p;
	int width, height;
public:
	EYXR(COE, src s);
	void cmp(COE, src, HX, int);
	void save2file();
};

class EYXU
{
public:
	float *p;
	int width, height;
public:
	EYXU(COE, src s);
	void cmp(COE, src, HX, int);
	void save2file();
};

class EYXD
{
public:
	float *p;
	int width, height;
public:
	EYXD(COE, src s);
	void cmp(COE, src, HX, int);
	void save2file();
};

#endif