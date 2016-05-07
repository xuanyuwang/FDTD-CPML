#include <iostream>
#include <string>
#include "E.h"
#include "H.h"
#include "cvl.h"
#include "source.h"
#include "computation.h"

#define CPML

using namespace std;

void main()
{
	src s(10, 10, 400);
	COE c(3, s);
	H h(c, s);
	E e(c, s);
	HXZ hxz(c, s);
	HYZ hyz(c, s);
	EXY exy(c, s);
	EYX eyx(c, s);

	int i, j, time;

	e.ez->checkout();
	//computation
	for (time = 0; time < s.size_time; time++)
	{
		cmp_hy(h, e, hxz, c, s);
		cmp_hx(h, e, hyz, c, s);
		cmp_ez(e, h, exy, eyx, c, s, time);
		//h.hy->save2file(h.hy->filename);
		//h.hx->save2file(h.hx->filename);
		s.cmp(time, &e.ez->p.at(e.ez->height / 2 * e.ez->width + e.ez->width / 2));
		e.ez->save2file(e.ez->filename);
	}
	//	for (j = 0; j < e.ez->width; j++)
	//	{
	//		s.cmp(i, &e.ez->p[e.ez->grid_num/2]);
	//	}
}