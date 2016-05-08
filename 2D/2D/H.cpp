#include "H.h"

H::H(COE c, src s)
{
	hy = new area(s.size_x + 2 * c.num_layer, "hy.txt");
	coe_h = s.dt / mu;
}

void H::real_pos_hy(int i, int j)
{
	cout << "Hy real position: (" << i+0.5 << ", " << j << ")" << endl;
}
