#pragma once
#ifndef Ez_H
#define Ez_H

#include "source.h"
#include "H.h"
#include "cvl.h"
#include <iostream>
#include <fstream>

using namespace std;

class E{
public:
	area *ez=NULL;
	float coe_E;
	//float coe_mur;
	const float epsilon = 8.85e-12f;
public:
	E(COE, src);
	void real_pos(int, int);	
};

#endif