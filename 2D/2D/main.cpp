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
	//Initialize the simulation area
	src s(10, 10, 300);
	COE coe(3, s);

	//creat an area for every part of field
	area *hx = new area(s.size_x + 2 * coe.num_layer, s.size_y + 2 * coe.num_layer + 1, "hx.txt");
	area *hy = new area(s.size_x + 2 * coe.num_layer + 1, s.size_y + 2 * coe.num_layer, "hy.txt");
	area *ez = new area(s.size_x + 2 * coe.num_layer + 1, s.size_y + 2 * coe.num_layer + 1, "ez.txt");
	area *hxz = new area(hy->width, hy->height, "hxz.txt");
	area *hyz = new area(hx->width, hx->height, "hyz.txt");
	area *exy = new area(ez->width, ez->height, "exy.txt");
	area *eyx = new area(ez->width, ez->height, "eyx.txt");

	//simulating
	int time;
	for (time = 0; time < s.size_time; time++)
	{
		//compute convolution terms Hyz & Hxz.
		cmp_hyz(hyz, ez, coe, s, time);
		cmp_hxz(hxz, ez, coe, s);
		
		//compute magnetics
		cmp_hx(hx, ez, hyz, coe, s);
		cmp_hy(hy, ez, hxz, coe, s);

		//compute convolution terms Exy & Eyx
		cmp_exy(exy, hy, coe, s);
		cmp_eyx(eyx, hx, coe, s, time);

		//compute electrics
		cmp_ez(ez, hy, hx, exy, eyx, coe, s, time);

		//assigned the value of source to eletrics
		s.cmp(time, &ez->p.at(ez->height / 2 * ez->width + ez->width / 2));

		//save data to files
		ez->save2file(ez->filename);
		//exy->save2file(exy->filename);
		//eyx->save2file(eyx->filename);
		//hyz->save2file(hyz->filename);
	}
	delete hx, hy, ez, exy, eyx, hyz, hxz;
}