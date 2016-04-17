#pragma once
#ifndef H_H
#define H_H

#include "source.h"
#include <iostream>
#include <fstream>
#include "E.h"
#include "cvl.h"

using namespace std;
class E;
class cvl;
class H{
public:
	float *Hy, coe_H, coe_H_cvl;
	int size_Hy;
	const float PI = 3.14159265f;
	const float mu = (4.0 * PI) * 1e-7f;

public:
	H(src, cvl);
	void checkout();
	void save2file();
	void cmp(E, cvl, src, int);
};

#endif