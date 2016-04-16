#pragma once
#ifndef EX_H
#define EX_H

#include "source.h"
#include "H.h"
#include <iostream>
#include <fstream>

using namespace std;

class H;
class cvl;
class E{
public:
	float *Ex, coe_E, coe_E_cvl, coe_mur;
	float Ex_nbdr, Ex_nbdl, Ex_bdr, Ex_bdl;
	int size_Ex;
	const float epsilon = 8.85e-12f;
	const float C = 3e8;
public:
	E(src source);
	void cmp(H Hy);
	void cmp(H Hy, cvl cvln);
	void checkout();
	void boundary();
	void save2file();
};

#endif