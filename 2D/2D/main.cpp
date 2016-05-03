#include <iostream>
#include <string>
#include "E.h"
#include "H.h"
#include "cvl.h"
#include "source.h"

#define CPML

using namespace std;

void main()
{
	src s(10, 10, 10);
	COE c(3, s);
	HXZL hxzl(c, s);
	HXZR hxzr(c, s);
	HXZU hxzu(c, s);
	HXZD hxzd(c, s);
	HYZL hyzl(c, s);
	HYZR hyzr(c, s);
	HYZU hyzu(c, s);
	HYZD hyzd(c, s);
	EXYL exyl(c, s);
	EXYR exyr(c, s);
	EXYU exyu(c, s);
	EXYD exyd(c, s);
	EYXL eyxl(c, s);
	EYXR eyxr(c, s);
	EYXU eyxu(c, s);
	EYXD eyxd(c, s);
	HX hx(c, s);
	HY hy(c, s);
	E e(s, c);

	for (int i = 0; i < s.size_time; i++)
	{
#ifdef CPML
		hxzl.cmp(c, s, e, i);
		hxzr.cmp(c, s, e, i);
		hxzu.cmp(c, s, e, i);
		hxzd.cmp(c, s, e, i);
		hyzl.cmp(c, s, e, i);
		hyzr.cmp(c, s, e, i);
		hyzu.cmp(c, s, e, i);
		hyzd.cmp(c, s, e, i);
		hy.cmp(e, c, s, hxzl, hxzr, hxzu, hxzd, i);
		hx.cmp(e, c, s, hyzl, hyzr, hyzu, hyzd, i);
		exyl.cmp(c, s, hy, i);
		exyr.cmp(c, s, hy, i);
		exyu.cmp(c, s, hy, i);
		exyd.cmp(c, s, hy, i);
		eyxl.cmp(c, s, hx, i);
		eyxr.cmp(c, s, hx, i);
		eyxu.cmp(c, s, hx, i);
		eyxd.cmp(c, s, hx, i);
		e.cmp(hx, hy, c, s,
			exyl, exyr, exyu, exyd,
			eyxl, eyxr, eyxu, eyxd,
			i);
#endif
		s.cmp(i, &e.Ez[e.num_grid / 2]);
		if (i < 10)
		{
			hy.save2file();
			hxzl.save2file();
			e.save2file();
			//cvln.save2file();
		}
	}
	//	//cvln.save2file_coe();
}