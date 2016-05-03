#include "H.h"
//#define debug

HX::HX(COE c, src s)
{
	width = s.size_x + 2 * c.num_layer;
	height = s.size_y + 2 * c.num_layer + 1;
	num_grid = width*height;
	coe_h = s.dt / (mu * s.dz);
	coe_h_cvl = s.dt / mu;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Hx.txt", ios::out);
	myfile.close();
}

void HX::cmp(E e, COE c, src s, HYZL hyzl, HYZR hyzr, HYZU hyzu, HYZD hyzd, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;

	//X: [0, 49], Y: [0, 50]
	int ysec1o = 0, ysec1c = pmlbd;//[0, 9]
	int ysec2o = pmlbd + 1, ysec2c = pmlbd + s.size_y + 1;//[10, 40]
	int ysec3o = height - c.num_layer, ysec3c = height - 1;//[41, 50]
	int xsec1o = 0, xsec1c = pmlbd;//[0,9]
	int xsec2o = pmlbd + 1, xsec2c = pmlbd + s.size_x;//[10, 39]
	int xsec3o = width - c.num_layer, xsec3c = width - 1;//[40, 49]

	//bottom area, include bottom-left CPML, bottom-mid CPML, and bottom-right CPML
	for (i = ysec1o; i <= ysec1c; i++)
	{
		//bottom-left CPML corner
		for (j = xsec1o; j <= xsec1c; j++)
		{
			p[i * width + j] += (-coe_h / c.kappa_half[pmlbd - j])
				*(e.Ez[i*e.size_x + j + 1] - e.Ez[i*e.size_x + j])
				+ coe_h_cvl * (-hyzl.p[i*hyzl.width + j]);
		}
		//bottom-mid CPML
		for (j = xsec2o; j <= xsec2c; j++)
		{
			p[i * width + j] += (-coe_h / c.kappa_full[pmlbd - i])
				*(e.Ez[i * e.size_x + j + 1] - e.Ez[i * e.size_x + j])
				+ coe_h_cvl*(-hyzd.p[i *hyzd.width + j - xsec2o]);
		}
		//bottom-right CPML
		for (j = xsec3o; j <= xsec3c; j++)
		{
			p[i * width + j] += (-coe_h / c.kappa_half[j - xsec3o])
				*(e.Ez[i * e.size_x + j + 1] - e.Ez[i * e.size_x + j])
				+ coe_h_cvl*(-hyzr.p[i*hyzr.width + j - xsec3o]);
		}
	}

	//mid area, include mid-left CPML, mid-right CPML, and simulation area
	for (i = ysec2o; i <= ysec2c; i++)
	{
		//mid-left CPML
		for (j = xsec1o; j <= xsec1c; j++)
		{
			p[i * width + j] += (-coe_h / c.kappa_half[pmlbd - j])
				*(e.Ez[i*e.size_x + j + 1] - e.Ez[i*e.size_x + j])
				+ coe_h_cvl * (-hyzl.p[i*hyzl.width + j]);
		}
		//simulation area
		for (j = xsec2o; j <= xsec2c; j++)
		{
			p[i * width + j] += -coe_h *(e.Ez[i*e.size_x + j + 1] - e.Ez[i*e.size_x + j]);
		}
		//mid-right CPML
		for (j = xsec3o; j <= xsec3c; j++)
		{
			p[i * width + j] += (-coe_h / c.kappa_half[j - xsec3o])
				*(e.Ez[i * e.size_x + j + 1] - e.Ez[i * e.size_x + j])
				+ coe_h_cvl*(-hyzr.p[i*hyzr.width + j - xsec3o]);
		}
	}

	//upper CPML.
	for (i = ysec3o; i <= ysec3c; i++)
	{
		//upper-left corner area
		for (j = xsec1o; j <= xsec1c; j++)
		{
			p[i * width + j] += (-coe_h / c.kappa_half[pmlbd - j])
				*(e.Ez[i*e.size_x + j + 1] - e.Ez[i*e.size_x + j])
				+ coe_h_cvl * (-hyzl.p[i*hyzl.width + j]);
		}
		//upper-center area
		for (j = xsec2o; j <= xsec2c; j++)
		{
			p[i * width + j] += (-coe_h / c.kappa_full[i - ysec3o])
				*(e.Ez[i * e.size_x + j + 1] - e.Ez[i * e.size_x + j])
				+ coe_h_cvl*(-hyzu.p[(i - ysec3o) *hyzu.width + j - xsec2o]);
		}
		//upper-right corner area
		for (j = xsec3o; j <= xsec3c; j++)
		{
			p[i * width + j] += (-coe_h / c.kappa_half[j - xsec3o])
				*(e.Ez[i * e.size_x + j + 1] - e.Ez[i * e.size_x + j])
				+ coe_h_cvl*(-hyzr.p[i*width + j - xsec3o]);
		}
	}
}

void HX::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Hx.txt", ios::app);
	for (int i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i * width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

HY::HY(COE c, src s)
{
	width = s.size_x + 2 * c.num_layer + 1;
	height = s.size_y + 2 * c.num_layer;
	num_grid = width*height;
	coe_h = s.dt / (mu * s.dz);
	coe_h_cvl = s.dt / mu;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Hy.txt", ios::out);
	myfile.close();
}

void HY::cmp(E e, COE c, src s, HXZL hxzl, HXZR hxzr, HXZU hxzu, HXZD hxzd, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;

	//X: [0, 50], Y: [0, 49]
	int ysec1o = 0, ysec1c = pmlbd;//[0, 9]
	int ysec2o = pmlbd + 1, ysec2c = pmlbd + s.size_y;//[10, 39]
	int ysec3o = height - c.num_layer, ysec3c = height - 1;//[40, 49]

	int xsec1o = 0, xsec1c = pmlbd;//[0,9]
	int xsec2o = pmlbd + 1, xsec2c = pmlbd + s.size_x + 1;//[10, 40]
	int xsec3o = width - c.num_layer, xsec3c = width - 1;//[41, 50]

//	//if (time == 0)
//	//{
//	//	cout << "computation of Hy" << endl;
//	//	cout << "y section 1:\t" << ysec1o << "\t" << ysec1c << endl;
//	//	cout << "y section 2:\t" << ysec2o << "\t" << ysec2c << endl;
//	//	cout << "y section 3:\t" << ysec3o << "\t" << ysec3c << endl;

//	//	cout << "x section 1:\t" << xsec1o << "\t" << xsec1c << endl;
//	//	cout << "x section 2:\t" << xsec2o << "\t" << xsec2c << endl;
//	//	cout << "x section 3:\t" << xsec3o << "\t" << xsec3c << endl;

//	//	cout << endl;
//	//}

	//bottom area, include bottom-left CPML, bottom-mid CPML, and bottom-right CPML
	for (i = ysec1o; i <= ysec1c; i++)
	{
		//bottom-left CPML corner
		for (j = xsec1o; j <= xsec1c; j++)
		{
			p[i * width + j] += (coe_h / c.kappa_full[pmlbd - j])
				*(e.Ez[(i + 1)*e.size_x + j] - e.Ez[i*e.size_x + j])
				+ coe_h_cvl * hxzl.p[i*hxzl.width + j];
		}
		//bottom-mid CPML
		for (j = xsec2o; j <= xsec2c; j++)
		{
			p[i * width + j] += (coe_h / c.kappa_half[pmlbd - i])
				*(e.Ez[(i + 1) * e.size_x + j] - e.Ez[i * e.size_x + j])
				+ coe_h_cvl*hxzd.p[i *hxzd.width + j - xsec2o];
		}
		//bottom-right CPML
		for (j = xsec3o; j <= xsec3c; j++)
		{
			p[i * width + j] += (coe_h / c.kappa_full[j - xsec3o])
				*(e.Ez[(i + 1) * e.size_x + j] - e.Ez[i * e.size_x + j])
				+ coe_h_cvl*hxzr.p[i*hxzr.width + j - xsec3o];
		}
	}

	//mid area, include mid-left CPML, mid-right CPML, and simulation area
	for (i = ysec2o; i <= ysec2c; i++)
	{
		//mid-left CPML
		for (j = xsec1o; j <= xsec1c; j++)
		{
			p[i * width + j] += (coe_h / c.kappa_full[pmlbd - j])
				*(e.Ez[(i + 1)*e.size_x + j] - e.Ez[i*e.size_x + j])
				+ coe_h_cvl * hxzl.p[i*c.side_sx + j];
		}
		//simulation area
		for (j = xsec2o; j <= xsec2c; j++)
		{
			p[i * width + j] += coe_h * (e.Ez[(i + 1)*e.size_x + j] - e.Ez[i*e.size_x + j]);
		}
		//mid-right CPML
		for (j = xsec3o; j <= xsec3c; j++)
		{
			p[i * width + j] += (coe_h / c.kappa_full[j - xsec3o])
				*(e.Ez[(i + 1) * e.size_x + j] - e.Ez[i * e.size_x + j])
				+ coe_h_cvl*hxzr.p[i*c.side_sx + j - xsec3o];
//			if (time == 7 && (i == 7 || i == 8) && j == 14)
//			{
//				cout << "Hy(" << i << ", " << j << "): " << Hy[i * size_Hy_x + j] << endl;
//				cout << "\tDebug\t"
//					<< "first term: " <<"Ez("<<i+1<<", "<<j+1<<")\t" <<Ez.Ez[(i + 1) * Ez.size_x + j] << ", \t"
//					<< //					Ez.Ez[i * Ez.size_x + j] << "\t"
//					//<< "second term: " << coe_//					H_cvl*c.Hxzr[i*c.side_sx + j - xsec3o] << "\t"
//					//					<< endl;
//					//			}
		}
	}

	//upper CPML.
	for (i = ysec3o; i <= ysec3c; i++)
	{
		//upper-left corner area
		for (j = xsec1o; j <= xsec1c; j++)
		{
			p[i * width + j] += (coe_h / c.kappa_full[j])*(e.Ez[(i + 1)*e.size_x + j] - e.Ez[i*e.size_x + j])
				+ coe_h_cvl * hxzl.p[i*hxzl.width + j];
		}
		//upper-center area
		for (j = xsec2o; j <= xsec2c; j++)
		{
			p[i * width + j] += (coe_h / c.kappa_half[i - ysec3o])*(e.Ez[(i + 1) * e.size_x + j] - e.Ez[i * e.size_x + j])
				+ coe_h_cvl*hxzu.p[(i - ysec3o) *hxzu.width + j - xsec2o];
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
			p[i * width + j] += (coe_h / c.kappa_full[j - xsec3o])*(e.Ez[(i + 1)* e.size_x + j] - e.Ez[i * e.size_x + j])
				+ coe_h_cvl*(-hxzr.p[i*hxzr.width + j - xsec3o]);
		}
	}
}

void HY::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Hy.txt", ios::app);
	for (int i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i * width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();


}

