#include <iostream>
#include <cmath>
#include <fstream>
#include "E.h"
#include "H.h"
#include "source.h"

using namespace std;

void main()
{
	src source(30, 300);
	H hy(source);
	E ex(source);

	source.checkout();
	hy.checkout();
	ex.checkout();

	for (int i = 0; i < source.size_time; i++)
	{
		hy.cmp(ex);
		ex.cmp(hy);
		ex.boundary();
		source.cmp(i, &ex.Ex[(int)(ex.size_Ex / 2)]);

		hy.save2file();
		ex.save2file();
	}
}