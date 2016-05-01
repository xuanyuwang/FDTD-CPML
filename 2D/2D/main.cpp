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
	src s(10, 10, 50);
	cvl cvln(3, s);
	H hy(s, cvln);
	E Ez(s, cvln);

	//s.checkout();
	//hy.checkout();
	//Ez.checkout();
	cvln.checkout();

	for (int i = 0; i < s.size_time; i++)
	{
#ifdef CPML
		cvln.cmp_cvlh(s, Ez, hy, i);
		hy.cmp_Hx(Ez, cvln, s, i);
		hy.cmp_Hy(Ez, cvln, s, i);
		cvln.cmp_cvle(s, Ez, hy, i);
		Ez.cmp(hy, cvln, s, i);
#endif
		s.cmp(i, &Ez.Ez[Ez.num_grid / 2]);
		if (i < 10)
		{
			hy.save2file();
			Ez.save2file();
			cvln.save2file();
		}
	}
	//cvln.save2file_coe();
}