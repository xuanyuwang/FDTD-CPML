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
	src s(30, 30, 300);
	cvl cvln(10, s);
	//H hy(s, cvln);
	//E ex(s, cvln);

	s.checkout();
	//hy.checkout();
	//ex.checkout();
	//cvln.checkout();

//	for (int i = 0; i < s.size_time; i++)
//	{
//#ifdef CPML
//		cvln.cmp_cvlh(s, ex, hy, i);
//		hy.cmp(ex, cvln, s, i);
//		cvln.cmp_cvle(s, ex, hy, i);
//		ex.cmp(hy, cvln, s, i);
//#else 
//		hy.cmp(ex);
//		ex.cmp(hy);
//		ex.boundary();
//#endif
//		s.cmp(i, &ex.Ex[ex.size_Ex / 2]);
//		hy.save2file();
//		ex.save2file();
//
//#ifdef CPML
//		cvln.save2file();
//#endif
//	}
	//cvln.save2file_coe();
}