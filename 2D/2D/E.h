#pragma once
#ifndef Ez_H
#define Ez_H

#include "source.h"
#include "H.h"
#include "cvl.h"
#include <iostream>
#include <fstream>

using namespace std;

class H;
class cvl;
class E{
public:
	float *Ez, coe_E, coe_E_cvl, coe_mur;
	int size_x, size_y;
	int num_grid;
	const float epsilon = 8.85e-12f;
	//const float C = 3e8;
public:
	E(src, cvl);
	void cmp(H Hy, cvl c, src s, int);
	void checkout();
	void save2file();
};

#endif