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
	src s(10, 10, 150);
	COE coe(3, s);
	H h(coe, s);
	E e(coe, s);
	HXZ hxz(coe, s);
	HYZ hyz(coe, s);
	EXY exy(coe, s);
	EYX eyx(coe, s);

	int i, j, time;

	//e.ez->checkout();
	//s.checkout();
	//coe.checkout(s);
	//computation
	for (time = 0; time < s.size_time; time++)
	{
		cmp_hyz(e, hyz, coe, s);
		cmp_hxz(e, hxz, coe, s);

		cmp_hy(h, e, hxz, coe, s);
		cmp_hx(h, e, hyz, coe, s);

		cmp_exy(h, exy, coe, s);
		cmp_eyx(h, eyx, coe, s);

		cmp_ez(e, h, exy, eyx, coe, s, time);
		//h.hy->save2file(h.hy->filename);
		//h.hx->save2file(h.hx->filename);
		s.cmp(time, &e.ez->p.at(e.ez->height / 2 * e.ez->width + e.ez->width / 2));
		//hyz.full->save2file(hyz.full->filename);
		e.ez->save2file(e.ez->filename);
	}
}