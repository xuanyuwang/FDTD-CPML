#pragma once
#ifndef Ez_H
#define Ez_H

#include "source.h"
#include "H.h"
#include "cvl.h"
#include <iostream>
#include <fstream>

using namespace std;

class HX;
class HY;
class COE;
class EXYL;
class EXYR;
class EXYU;
class EXYD;
class EYXL;
class EYXR;
class EYXU;
class EYXD;

class E{
public:
	float *Ez, coe_E, coe_E_cvl, coe_mur;
	int size_x, size_y;
	int num_grid;
	const float epsilon = 8.85e-12f;
	//const float C = 3e8;
public:
	E(src, COE);
	void E::cmp(HX, HY, COE c, src s,
		EXYL, EXYR, EXYU , EXYD ,
		EYXL , EYXR, EYXU , EYXD ,
		int);
	void checkout();
	void save2file();
};

#endif