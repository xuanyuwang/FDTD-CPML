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
	src source(30, 100);
	H hy(source);
	E ex(source);
	cvl cvln(10, source);

	//source.checkout();
	//hy.checkout();
	//ex.checkout();
	//cvln.checkout();

	for (int i = 0; i < source.size_time; i++)
	{
#ifdef CPML
		hy.cmp(ex);
		ex.cmp(hy,cvln);
		cvln.cmp_H_l(source, ex.Ex[0]);
		cvln.cmp_H_r(source, ex.Ex[ex.size_Ex-1]);
		cvln.cmp_Hy(source, ex.Ex[0], ex.Ex[ex.size_Ex-1]);
		cvln.cmp_E_l(source);
		cvln.cmp_E_r(source);
		cvln.cmp_Ex(source);
#else 
		hy.cmp(ex);
		ex.cmp(hy);
		ex.boundary();
#endif
		source.cmp(i, &ex.Ex[ex.size_Ex / 2]);
		hy.save2file();
		ex.save2file();

#ifdef CPML
		cvln.save2file();
#endif
	}
	cvln.save2file_coe();
}