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
	src s(10, 150);
	COE coe(3, s);
	H h(coe, s);
	E e(coe, s);
	area *hzx = new area(h.hy->length, "hzx.txt");
	area *ezy = new area(e.ex->length, "ezy.txt");
	int i, j, time;

	//e.ez->checkout();
	//s.checkout();
	coe.checkout(s);
	//computation
	for (time = 0; time < s.size_time; time++)
	{
		cmp_hzx(e, hzx, coe, s, time);
		//hzx->save2file(hzx->filename);
		cmp_hy(h, e, hzx, coe, s, time);
		//h.hy->save2file(h.hy->filename);
		cmp_ezy(h, ezy, coe, s);
		//ezy->save2file(ezy->filename);
		cmp_ex(e, h, ezy, coe, s, time);
		s.cmp(time, &e.ex->p.at(e.ex->length / 2));
		e.ex->save2file(e.ex->filename);
	}
}