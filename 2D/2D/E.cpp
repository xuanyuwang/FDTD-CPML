#include "E.h"
//#define MUR
//#define debug

E::E(COE c,src s)
{
	ex = new area(s.size_x + 2 * c.num_layer + 1, "ex.txt");
	coe_E = s.dt / epsilon;
	//coe_mur = (C*s.dt - s.dz) / (C*s.dt + s.dz);
}

void E::real_pos(int i, int j)
{
	cout << "real position: (" << i << ", " << j << ")" << endl;
}