#include "H.h"
//#define debug

H::H(COE c, src s)
{
	hx = new area(s.size_x + 2 * c.num_layer, s.size_y + 2 * c.num_layer + 1, "hx.txt");
	hy = new area(s.size_x + 2 * c.num_layer + 1, s.size_y + 2 * c.num_layer, "hy.txt");
	coe_h = s.dt / mu;
}

void H::real_pos_hx(int i, int j)
{
	cout << "Hx real position: (" << i << ", " << j + 0.5 << ")" << endl;
}

void H::real_pos_hy(int i, int j)
{
	cout << "Hy real position: (" << i+0.5 << ", " << j << ")" << endl;
}
