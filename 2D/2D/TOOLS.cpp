#include "TOOLS.h"

void TOOLS::alloc(float* &p, int x, int y)
{
	p = (float*)malloc(x*y*sizeof(float));
	for (int i=0; i < x*y;i++)
	{
		p[i] = 0.f;
	}
}

TOOLS::TOOLS()
{
}


TOOLS::~TOOLS()
{
}
