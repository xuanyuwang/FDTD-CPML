#pragma once
#ifndef EX_H
#define EX_H

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
	float *Ex, coe_E, coe_E_cvl, coe_mur;
	int size_Ex;
	const float epsilon = 8.85e-12f;
	//const float C = 3e8;
public:
	E(src, cvl);
	void cmp(H Hy, cvl c, src s, int);
	void checkout();
	void save2file();
};

#endif