#pragma once
#ifndef H_H
#define H_H

#include "source.h"
#include <iostream>
#include <fstream>
#include "E.h"
#include "cvl.h"
#include "TOOLS.h"

using namespace std;
class E;
class HYZL;
class HXZL;
class HXZR;
class HXZU;
class HXZD;
class HYZR;
class HYZU;
class HYZD;
class COE;

class HX
{
public:
	int width, height;
	int num_grid;
	float *p;
	float coe_h, coe_h_cvl;
	const float PI = 3.14159265f;
	const float mu = (4.0 * PI) * 1e-7f;

public:
	HX(COE, src);
	void cmp(E, COE, src, HYZL, HYZR, HYZU, HYZD, int);
	void save2file();
};

class HY
{
public:
	int width, height;
	int num_grid;
	float *p;
	float coe_h, coe_h_cvl;
	const float PI = 3.14159265f;
	const float mu = (4.0 * PI) * 1e-7f;

public:
	HY(COE, src);
	void cmp(E e, COE c, src s, HXZL hxzl, HXZR hxzr, HXZU hxzu, HXZD hxzd, int time);
	void save2file();
};

#endif