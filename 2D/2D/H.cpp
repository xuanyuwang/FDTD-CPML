#include "H.h"
//#define debug

H::H(src s, cvl c)
{
	int i;

	size_Hx_x = s.size_x + 2 * c.num_layer;
	size_Hx_y = s.size_y + 2 * c.num_layer + 1;
	num_grid_Hx = size_Hx_y * size_Hx_x;

	size_Hy_x = s.size_x + 2 * c.num_layer + 1;
	size_Hy_y = s.size_y + 2 * c.num_layer;
	num_grid_Hy = size_Hy_y * size_Hy_x;

	Hy = (float *)malloc(num_grid_Hy * sizeof(float));
	memset(Hy, 0, num_grid_Hy * sizeof(float));
	Hx = (float *)malloc(num_grid_Hx * sizeof(float));
	memset(Hx, 0, num_grid_Hx * sizeof(float));

	coe_H = s.dt / (mu * s.dz);
	coe_H_cvl = s.dt / mu;

	fstream myfile;
	myfile.open("Hy.txt", ios::out);
	myfile.close();

	myfile.open("Hx.txt", ios::out);
	myfile.close();
}

void H::cmp_Hx(E Ez, cvl c, src s, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;

	//X: [0, 49], Y: [0, 50]
	int ysec1o = 0, ysec1c = pmlbd;//[0, 9]
	int ysec2o = pmlbd + 1, ysec2c = pmlbd + s.size_y + 1;//[10, 40]
	int ysec3o = size_Hx_y - c.num_layer, ysec3c = size_Hx_y - 1;//[41, 50]
	int xsec1o = 0, xsec1c = pmlbd;//[0,9]
	int xsec2o = pmlbd + 1, xsec2c = pmlbd + s.size_x;//[10, 39]
	int xsec3o = size_Hx_x - c.num_layer, xsec3c = size_Hx_x - 1;//[40, 49]

	//bottom area, include bottom-left CPML, bottom-mid CPML, and bottom-right CPML
	for (i = ysec1o; i <= ysec1c; i++)
	{
		//bottom-left CPML corner
		for (j = xsec1o; j <= xsec1c; j++)
		{
			Hx[i * size_Hx_x + j] += (-coe_H / c.kappa_half[pmlbd - j])*(Ez.Ez[i*Ez.size_x + j + 1] - Ez.Ez[i*Ez.size_x + j])
				+ coe_H_cvl * (-c.Hyzl[i*c.side_sx + j]);
		}
		//bottom-mid CPML
		for (j = xsec2o; j <= xsec2c; j++)
		{
			Hx[i * size_Hx_x + j] += (-coe_H / c.kappa_full[pmlbd - i])*(Ez.Ez[i * Ez.size_x + j + 1] - Ez.Ez[i * Ez.size_x + j])
				+ coe_H_cvl*(-c.Hyzd[i *c.side_sx + j - xsec2o]);
		}
		//bottom-right CPML
		for (j = xsec3o; j <= xsec3c; j++)
		{
			Hx[i * size_Hx_x + j] += (-coe_H / c.kappa_half[j - xsec3o])*(Ez.Ez[i * Ez.size_x + j + 1] - Ez.Ez[i * Ez.size_x + j])
				+ coe_H_cvl*(-c.Hyzr[i*c.side_sx + j - xsec3o]);
		}
	}

	//mid area, include mid-left CPML, mid-right CPML, and simulation area
	for (i = ysec2o; i <= ysec2c; i++)
	{
		//mid-left CPML
		for (j = xsec1o; j <= xsec1c; j++)
		{
			Hx[i * size_Hx_x + j] += (-coe_H / c.kappa_half[pmlbd - j])*(Ez.Ez[i*Ez.size_x + j + 1] - Ez.Ez[i*Ez.size_x + j])
				+ coe_H_cvl * (-c.Hyzl[i*c.side_sx + j]);
		}
		//simulation area
		for (j = xsec2o; j <= xsec2c; j++)
		{
			Hx[i * size_Hx_x + j] += -coe_H * (Ez.Ez[i*Ez.size_x + j + 1] - Ez.Ez[i*Ez.size_x + j]);
//			if (time == 5 && (i == 8) && j == 12)
//			{
//				cout << "Hx(" << i << ", " << j << "): " << Hx[i * size_Hx_x + j] << endl;
//				cout << "\tDebug\t"
//					<< "(" << i << ", " << j + 1 << ")\t"
//					<< Ez.Ez[i*Ez.size_x + j + 1]<<"\t"
//					<< "(" << i << ", " << j << ")\t"
//					<< Ez.Ez[i*Ez.size_x + j]
//					<< endl;
//			}
		}
		//mid-right CPML
		for (j = xsec3o; j <= xsec3c; j++)
		{
			Hx[i * size_Hx_x + j] += (-coe_H / c.kappa_half[j - xsec3o])*(Ez.Ez[i * Ez.size_x + j + 1] - Ez.Ez[i * Ez.size_x + j])
				+ coe_H_cvl*(-c.Hyzr[i*c.side_sx + j - xsec3o]);
		}
	}

	//upper CPML.
	for (i = ysec3o; i <= ysec3c; i++)
	{
		//upper-left corner area
		for (j = xsec1o; j <= xsec1c; j++)
		{
			Hx[i * size_Hx_x + j] += (-coe_H / c.kappa_half[pmlbd - j])*(Ez.Ez[i*Ez.size_x + j + 1] - Ez.Ez[i*Ez.size_x + j])
				+ coe_H_cvl * (-c.Hyzl[i*c.side_sx + j]);
		}
		//upper-center area
		for (j = xsec2o; j <= xsec2c; j++)
		{
			Hx[i * size_Hx_x + j] += (-coe_H / c.kappa_full[i - ysec3o])*(Ez.Ez[i * Ez.size_x + j + 1] - Ez.Ez[i * Ez.size_x + j])
				+ coe_H_cvl*(-c.Hyzu[(i - ysec3o) *c.side_sx + j - xsec2o]);
		}
		//upper-right corner area
		for (j = xsec3o; j <= xsec3c; j++)
		{
			Hx[i * size_Hx_x + j] += (-coe_H / c.kappa_half[j - xsec3o])*(Ez.Ez[i * Ez.size_x + j + 1] - Ez.Ez[i * Ez.size_x + j])
				+ coe_H_cvl*(-c.Hyzr[i*c.side_sx + j - xsec3o]);
		}
	}
}

void H::cmp_Hy(E Ez, cvl c, src s, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;

	//X: [0, 50], Y: [0, 49]
	int ysec1o = 0, ysec1c = pmlbd;//[0, 9]
	int ysec2o = pmlbd + 1, ysec2c = pmlbd + s.size_y;//[10, 39]
	int ysec3o = size_Hy_y - c.num_layer, ysec3c = size_Hy_y - 1;//[40, 49]

	int xsec1o = 0, xsec1c = pmlbd;//[0,9]
	int xsec2o = pmlbd + 1, xsec2c = pmlbd + s.size_x + 1;//[10, 40]
	int xsec3o = size_Hy_x - c.num_layer, xsec3c = size_Hy_x - 1;//[41, 50]

	//if (time == 0)
	//{
	//	cout << "computation of Hy" << endl;
	//	cout << "y section 1:\t" << ysec1o << "\t" << ysec1c << endl;
	//	cout << "y section 2:\t" << ysec2o << "\t" << ysec2c << endl;
	//	cout << "y section 3:\t" << ysec3o << "\t" << ysec3c << endl;

	//	cout << "x section 1:\t" << xsec1o << "\t" << xsec1c << endl;
	//	cout << "x section 2:\t" << xsec2o << "\t" << xsec2c << endl;
	//	cout << "x section 3:\t" << xsec3o << "\t" << xsec3c << endl;

	//	cout << endl;
	//}

	//bottom area, include bottom-left CPML, bottom-mid CPML, and bottom-right CPML
	for (i = ysec1o; i <= ysec1c; i++)
	{
		//bottom-left CPML corner
		for (j = xsec1o; j <= xsec1c; j++)
		{
			Hy[i * size_Hy_x + j] += (coe_H / c.kappa_full[pmlbd - j])*(Ez.Ez[(i + 1)*Ez.size_x + j] - Ez.Ez[i*Ez.size_x + j])
				+ coe_H_cvl * c.Hxzl[i*c.side_sx + j];
		}
		//bottom-mid CPML
		for (j = xsec2o; j <= xsec2c; j++)
		{
			Hy[i * size_Hy_x + j] += (coe_H / c.kappa_half[pmlbd - i])*(Ez.Ez[(i + 1) * Ez.size_x + j] - Ez.Ez[i * Ez.size_x + j])
				+ coe_H_cvl*c.Hxzd[i *(c.ud_sx + 1) + j - xsec2o];
		}
		//bottom-right CPML
		for (j = xsec3o; j <= xsec3c; j++)
		{
			Hy[i * size_Hy_x + j] += (coe_H / c.kappa_full[j - xsec3o])*(Ez.Ez[(i + 1) * Ez.size_x + j] - Ez.Ez[i * Ez.size_x + j])
				+ coe_H_cvl*c.Hxzr[i*c.side_sx + j - xsec3o];
		}
	}

	//mid area, include mid-left CPML, mid-right CPML, and simulation area
	for (i = ysec2o; i <= ysec2c; i++)
	{
		//mid-left CPML
		for (j = xsec1o; j <= xsec1c; j++)
		{
			Hy[i * size_Hy_x + j] += (coe_H / c.kappa_full[pmlbd - j])*(Ez.Ez[(i + 1)*Ez.size_x + j] - Ez.Ez[i*Ez.size_x + j])
				+ coe_H_cvl * c.Hxzl[i*c.side_sx + j];
		}
		//simulation area
		for (j = xsec2o; j <= xsec2c; j++)
		{
			Hy[i * size_Hy_x + j] += coe_H * (Ez.Ez[(i + 1)*Ez.size_x + j] - Ez.Ez[i*Ez.size_x + j]);
		}
		//mid-right CPML
		for (j = xsec3o; j <= xsec3c; j++)
		{
			Hy[i * size_Hy_x + j] += (coe_H / c.kappa_full[j - xsec3o])*(Ez.Ez[(i + 1) * Ez.size_x + j + 1] - Ez.Ez[i * Ez.size_x + j])
				+ coe_H_cvl*c.Hxzr[i*c.side_sx + j - xsec3o];
		}
	}

	//upper CPML.
	for (i = ysec3o; i <= ysec3c; i++)
	{
		//upper-left corner area
		for (j = xsec1o; j <= xsec1c; j++)
		{
			Hy[i * size_Hy_x + j] += (coe_H / c.kappa_full[j])*(Ez.Ez[(i + 1)*Ez.size_x + j] - Ez.Ez[i*Ez.size_x + j])
				+ coe_H_cvl * c.Hxzl[i*c.side_sx + j];
		}
		//upper-center area
		for (j = xsec2o; j <= xsec2c; j++)
		{
			Hy[i * size_Hy_x + j] += (coe_H / c.kappa_half[i - ysec3o])*(Ez.Ez[(i + 1) * Ez.size_x + j] - Ez.Ez[i * Ez.size_x + j])
				+ coe_H_cvl*c.Hxzu[(i - ysec3o) *(c.ud_sx + 1) + j - xsec2o];
			//if (time == 0 && i == ysec3c && j == 11)
			//{
			//	cout << "coordinates\t" << i << "\t" << j << endl;
			//	cout << "\t" << Hy[i*size_Hy_x + j] << endl;
			//	cout << "\t" << (i - ysec3o) *(c.ud_sx + 1) + j - xsec2o << "\t" << c.Hxzu[(i - ysec3o) *(c.ud_sx + 1) + j - xsec2o] << endl;
			//}
		}
		//upper-right corner area
		for (j = xsec3o; j <= xsec3c; j++)
		{
			Hy[i * size_Hy_x + j] += (coe_H / c.kappa_full[j - xsec3o])*(Ez.Ez[(i + 1)* Ez.size_x + j] - Ez.Ez[i * Ez.size_x + j])
				+ coe_H_cvl*(-c.Hxzr[i*c.side_sx + j - xsec3o]);
		}
	}
}

void H::checkout()
{
	//cout << "Hy: " << endl;
	//for (int i = 0; i < size_Hy; i++){
	//	cout << Hy[i] << "\t";
	//}
	//cout << endl;
}

void H::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Hy.txt", ios::app);
	for (int i = size_Hy_y - 1; i >= 0; i--)
	{
		for (j = 0; j < size_Hy_x; j++)
		{
			myfile << Hy[i * size_Hy_x + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();

	myfile.open("Hx.txt", ios::app);
	for (int i = size_Hx_y - 1; i >= 0; i--)
	{
		for (j = 0; j < size_Hx_x; j++)
		{
			myfile << Hx[i * size_Hx_x + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}