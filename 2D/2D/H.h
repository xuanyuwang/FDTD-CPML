#pragma once

#include "source.h"
#include <iostream>
#include <fstream>
#include "area.h"
#include "cvl.h"

using namespace std;

class H
{
public:
	area *hx=NULL, *hy=NULL;
	float coe_h;
	const float PI = 3.14159265f;
	const float mu = (4.0 * PI) * 1e-7f;
public:
	H(COE, src);
	void real_pos_hx(int, int);
	void real_pos_hy(int, int);
};