#include "E.h"
#define MUR
//#define debug

E::E(src s, COE c)
{
	size_x = s.size_x + 1 + 2 * c.num_layer;
	size_y = s.size_y + 1 + 2 * c.num_layer;
	num_grid = size_x * size_y;
	Ez = (float *)malloc(num_grid * sizeof(float));
	memset(Ez, 0, num_grid * sizeof(float));

	coe_E = s.dt / (epsilon * s.dz);
	coe_E_cvl = s.dt / epsilon;
	//coe_mur = (C*s.dt - s.dz) / (C*s.dt + s.dz);

	fstream myfile;
	myfile.open("Ez.txt", ios::out);
	myfile.close();
}

void E::cmp(HX hx, HY hy, COE c, src s,
	EXYL exyl, EXYR exyr, EXYU exyu, EXYD exyd,
	EYXL eyxl, EYXR eyxr, EYXU eyxu, EYXD eyxd,
	int time)
{
	int i, j;
	//[0,pmlbd]([0,9]), [pmlbd,pmlbd+s.size_spzce]([10, 40]), [size_Ez - 10, size_Ez - 1]([41,50])
	int pmlbd = c.num_layer - 1;
	int ysec1o = 0, ysec1c = pmlbd;
	int ysec2o = pmlbd + 1, ysec2c = pmlbd + s.size_y + 1;
	int ysec3o = size_y - c.num_layer, ysec3c = size_y - 1;

	int xsec1o = 0, xsec1c = pmlbd;
	int xsec2o = pmlbd + 1, xsec2c = pmlbd + s.size_x + 1;
	int xsec3o = size_x - c.num_layer, xsec3c = size_x - 1;
	int width = c.ud_sx + 1;

	//for bottom CPML
	for (i = ysec1o; i <= ysec1c; i++)
	{
		//left bottom corner, CPML
		for (j = xsec1o; j <= xsec1c; j++)
		{
			if (i == 0 || j == 0)
			{
				Ez[i*size_x + j] = 0.f;
				continue;
			}
			Ez[i*size_x + j] += (coe_E / c.kappa_full[pmlbd - j])*(
				(hy.p[i*hy.width + j] - hy.p[(i - 1)*hy.width + j]) -
				(hx.p[i*hx.width + j] - hx.p[i*hx.width + j - 1])
				)
				+ coe_E_cvl*(exyl.p[i*c.side_sx + j] - eyxl.p[i*c.side_sx + j]);
		}

		//center bottom corner, CPML
		for (j = xsec2o; j <= xsec2c; j++)
		{
			if (i == 0)
			{
				Ez[i*size_x + j] = 0.f;
				continue;
			}
			Ez[i*size_x + j] += (coe_E / c.kappa_full[pmlbd - i])*(
				(hy.p[i*hy.width + j] - hy.p[(i - 1)*hy.width + j]) -
				(hx.p[i*hx.width + j] - hx.p[i*hx.width + j - 1])
				)
				+ coe_E_cvl*(exyd.p[i*width + j - xsec2o] - eyxd.p[i*width + j - xsec2o]);
		}

		//right bottom corner, CPML
		for (j = xsec3o; j <= xsec3c; j++)
		{
			if (i == 0 || j == xsec3c)
			{
				Ez[i*size_x + j] = 0.f;
				continue;
			}
			Ez[i*size_x + j] += (coe_E / c.kappa_full[j - xsec3o])*(
				(hy.p[i*hy.width + j] - hy.p[(i - 1)*hy.width + j]) -
				(hx.p[i*hx.width + j] - hx.p[i*hx.width + j - 1])
				)
				+ coe_E_cvl*(exyr.p[i *c.side_sx + j - xsec3o] - eyxr.p[i*c.side_sx + j - xsec3o]);
		}
	}

	//mid
	for (i = ysec2o; i <= ysec2c; i++)
	{
		//mid left CPML
		for (j = xsec1o; j <= xsec1c; j++)
		{
			if (j == 0)
			{
				Ez[i*size_x + j] = 0.f;
				continue;
			}
			Ez[i*size_x + j] += (coe_E / c.kappa_full[pmlbd - j])*(
				(hy.p[i*hy.width + j] - hy.p[(i - 1)*hy.width + j]) -
				(hx.p[i*hx.width + j] - hx.p[i*hx.width + j - 1])
				)
				+ coe_E_cvl*(exyl.p[i*c.side_sx + j] - eyxl.p[i*c.side_sx + j]);
		}

		//mid center FDTD area
		for (j = xsec2o; j <= xsec2c; j++)
		{
			Ez[i*size_x + j] += coe_E*(
				(hy.p[i*hy.width + j] - hy.p[(i - 1)*hy.width + j]) -
				(hx.p[i*hx.width + j] - hx.p[i*hx.width + j - 1])
				);
		}

		//mid right CPML
		for (j = xsec3o; j <= xsec3c; j++)
		{
			if (j == xsec3c)
			{
				Ez[i*size_x + j] = 0.f;
				continue;
			}
			Ez[i*size_x + j] += (coe_E / c.kappa_full[j - xsec3o])*(
				(hy.p[i*hy.width + j] - hy.p[(i - 1)*hy.width + j]) -
				(hx.p[i*hx.width + j] - hx.p[i*hx.width + j - 1])
				)
				+ coe_E_cvl*(exyr.p[i*c.side_sx + j - xsec3o] - eyxr.p[i*c.side_sx + j - xsec3o]);
		}
	}

	//upper area
	for (i = ysec3o; i <= ysec3c; i++)
	{
		//upper left CPML corner
		for (j = xsec1o; j <= xsec1c; j++)
		{
			if (i == ysec3c || j == 0)
			{
				Ez[i*size_x + j] = 0.f;
				continue;
			}
			Ez[i*size_x + j] += (coe_E / c.kappa_full[pmlbd - j])*(
				(hy.p[i*hy.width + j] - hy.p[(i - 1)*hy.width + j]) -
				(hx.p[i*hx.width + j] - hx.p[i*hx.width + j - 1])
				)
				+ coe_E_cvl*(exyl.p[i*c.side_sx + j] - eyxl.p[i*c.side_sx + j]);
		}

		//upper mid CPML corner
		for (j = xsec2o; j <= xsec2c; j++)
		{
			if (i == ysec3c)
			{
				Ez[i*size_x + j] = 0.f;
				continue;
			}
			Ez[i*size_x + j] += (coe_E / c.kappa_full[i - ysec3o])*(
				(hy.p[i*hy.width + j] - hy.p[(i - 1)*hy.width + j]) -
				(hx.p[i*hx.width + j] - hx.p[i*hx.width + j - 1])
				)
				+ coe_E_cvl*(exyu.p[(i - ysec3o)*width + j - xsec2o] - eyxu.p[(i - ysec3o)*width + j - xsec2o]);
		}

		//upper right CPML corner
		for (j = xsec3o; j <= xsec3c; j++)
		{
			if (i == ysec3c || j == xsec3c)
			{
				Ez[i*size_x + j] = 0.f;
				continue;
			}
			Ez[i*size_x + j] += (coe_E / c.kappa_full[j - xsec3o])*(
				(hy.p[i*hy.width + j] - hy.p[(i - 1)*hy.width + j]) -
				(hx.p[i*hx.width + j] - hx.p[i*hx.width + j - 1])
				)
				+ coe_E_cvl*(exyr.p[i*c.side_sx + j - xsec3o] - eyxr.p[i*c.side_sx + j - xsec3o]);

		}
	}
}

void E::checkout()
{
}

void E::save2file()
{
	int i, j;
	fstream myfile;
	myfile.open("Ez.txt", ios::app);

	for (i = size_y - 1; i >= 0; i--)
	{
		for (j = 0; j < size_x; j++)
		{
			myfile << Ez[i*size_x + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}