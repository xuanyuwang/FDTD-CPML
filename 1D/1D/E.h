#pragma once
#ifndef EX_H
#define EX_H

#include "source.h"
#include "H.h"
#include <iostream>
#include <fstream>

using namespace std;

class H;
class E{
public:
	float *Ex, coe_E;
	int size_Ex;
	const float epsilon = 8.85e-12f;
public:
	E(src source);
	void cmp(H Hy);
	void checkout();
	void boundary();
	void save2file();
};

#endif