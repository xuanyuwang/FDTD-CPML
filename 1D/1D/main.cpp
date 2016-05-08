#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "E.h"
#include "H.h"
#include "cvl.h"
#include "source.h"

#define CPML

using namespace std;

void main()
{
	src s(10, 150);
	cvl cvln(3, s);
	H hy(s, cvln);
	E ex(s, cvln);

	//s.checkout();
	//hy.checkout();
	//ex.checkout();
	//cvln.checkout();

	for (int i = 0; i < s.size_time; i++)
	{
#ifdef CPML
		cvln.cmp_cvlh(s, ex, hy, i);
		hy.cmp(ex, cvln, s, i);
		//hy.save2file();
		cvln.cmp_cvle(s, ex, hy, i);
		//cvln.save2file();
		ex.cmp(hy, cvln, s, i);
#else 
		hy.cmp(ex);
		ex.cmp(hy);
		ex.boundary();
#endif
		s.cmp(i, &ex.Ex[ex.size_Ex / 2]);
		//ex.save2file(s,cvln);

#ifdef CPML
#endif
	}
	cvln.save2file_coe();
}