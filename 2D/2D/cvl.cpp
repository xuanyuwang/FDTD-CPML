#include "cvl.h"
#define debug

cvl::cvl(int number, src s)
{
	num_layer = number;
	d = num_layer * s.dz;
	sigma_max = (m + 1) / (sqrt(eps_r) * 150 * PI * s.dz);
	side_sx = num_layer;
	side_sy = s.size_y + 2 * num_layer;
	ud_sx = s.size_x;//10
	ud_sy = num_layer;//3

	alloc(Hxzl, side_sx, side_sy);//3, 16, for Hy
	alloc(Hyzl, side_sx, side_sy + 1);//3, 17, for Hx
	alloc(Exyl, side_sx, side_sy + 1);//3, 17
	alloc(Eyxl, side_sx, side_sy + 1);//3, 17

	alloc(Hxzr, side_sx, side_sy);
	alloc(Hyzr, side_sx, side_sy + 1);
	alloc(Exyr, side_sx, side_sy + 1);
	alloc(Eyxr, side_sx, side_sy + 1);//3,17

	alloc(Hxzu, ud_sx + 1, ud_sy);//11, 3
	alloc(Hyzu, ud_sx, ud_sy);//10, 3
	alloc(Exyu, ud_sx + 1, ud_sy);//11, 3
	alloc(Eyxu, ud_sx + 1, ud_sy);

	alloc(Hxzd, ud_sx + 1, ud_sy);
	alloc(Hyzd, ud_sx, ud_sy);
	alloc(Exyd, ud_sx + 1, ud_sy);
	alloc(Eyxd, ud_sx + 1, ud_sy);

	set_distance(s);
	set_kappa(s);
	set_sigma(s);
	set_alpha(s);
	set_c(s);
	set_full_coe(s);
	set_half_coe(s);

	int i;
	fstream myfile;
	//myfile.open("Hzx.txt", ios::out);
	//myfile.close();

	myfile.open("Exyl.txt", ios::out);
	myfile.close();

	myfile.open("Exyr.txt", ios::out);
	myfile.close();

	myfile.open("Eyxl.txt", ios::out);
	myfile.close();

	myfile.open("Eyxr.txt", ios::out);
	myfile.close();

	myfile.open("Hxzl.txt", ios::out);
	myfile.close();

	myfile.open("Hxzu.txt", ios::out);
	myfile.close();

	myfile.open("Hyzl.txt", ios::out);
	myfile.close();

	myfile.open("Hyzu.txt", ios::out);
	myfile.close();

	myfile.open("Hyzr.txt", ios::out);
	myfile.close();

	myfile.open("coe.txt", ios::out);
	myfile.close();
}


cvl::~cvl()
{
	//free(Hxzl);
	//free(Hxzr);
	//free(Hxzu);
	//free(Hxzd);

	//free(Hyzl);
	//free(Hyzr);
	//free(Hyzu);
	//free(Hyzd);

	//free(Ezyl);
	//free(Ezyr);
	//free(Ezyu);
	//free(Ezyd);

	//free(Eyxl);
	//free(Eyxr);
	//free(Eyxu);
	//free(Eyxd);

	//free(distance_full);
	//free(distance_half);
	//free(sigma_full);
	//free(sigma_half);
	//free(alpha_full);
	//free(alpha_half);
	////free(kappa_full);
	//free(kappa_half);
	//free(c_E);
	//free(c_H);
	//free(cvl_H_coe);
	//free(cvl_E_coe);
}

void cvl::alloc(float* &p, int x, int y)
{
	p = (float*)malloc(x * y * sizeof(float));
	for (int i = 0; i < x*y; i++)
	{
		p[i] = 0.f;
	}
	//memset(p, 0, x*y*sizeof(float));
}

//the distance of every field points of H and E in CPML layer.
//the "i_min" is the nearest distance, the "I_max" is the furthest distance.
void cvl::set_distance(src source)
{
	int i;
	distance_full = (float *)malloc(num_layer*sizeof(float));
	distance_half = (float *)malloc(num_layer*sizeof(float));
	//for distance_full
	for (i = 0; i < num_layer; i++){
		distance_full[i] = (i + 1) * source.dz;
	}
	//for distance_half
	for (i = 0; i < num_layer; i++){
		distance_half[i] = (i + 0.5) * source.dz;
	}
}

void cvl::set_sigma(src source)
{
	int i;
	sigma_full = (float *)malloc(num_layer * sizeof(float));
	sigma_half = (float *)malloc(num_layer * sizeof(float));
	//sigma E
	for (i = 0; i < num_layer; i++)
	{
		sigma_full[i] = sigma_max * powf((distance_full[i] / d), m);
	}

	//sigma_half
	for (i = 0; i < num_layer; i++)
	{
		sigma_half[i] = sigma_max * powf((distance_half[i] / d), m);
	}
}

void cvl::set_kappa(src source)
{
	kappa_full = (float *)malloc(num_layer * sizeof(float));
	kappa_half = (float *)malloc(num_layer * sizeof(float));
	int i;
	//for kappa_full
	for (i = 0; i < num_layer; i++)
	{
		kappa_full[i] = 1 + (kappa_max - 1) * powf((distance_full[i] / d), m);
	}
	//for kappa_half
	for (i = 0; i < num_layer; i++)
	{
		kappa_half[i] = 1 + (kappa_max - 1) * powf((distance_half[i] / d), m);
	}
}

void cvl::set_alpha(src source)
{
	alpha_full = (float *)malloc(num_layer * sizeof(float));
	alpha_half = (float *)malloc(num_layer * sizeof(float));

	int i;
	//for alpha_full
	for (i = 0; i < num_layer; i++)
	{
		alpha_full[i] = sigma_full[i] / (eps0 * kappa_full[i]);
	}
	//for alpha_half
	for (i = 0; i < num_layer; i++)
	{
		alpha_half[i] = sigma_half[i] / (eps0 * kappa_half[i]);
	}
}

void cvl::set_c(src source)
{
	c_full = (float *)malloc(num_layer * sizeof(float));
	c_half = (float *)malloc(num_layer * sizeof(float));

	int i;
	//for E
	for (i = 0; i < num_layer; i++)
	{
		c_full[i] = (1 / kappa_full[i])*(exp(-alpha_full[i] * source.dt) - 1);
	}
	//for H
	for (i = 0; i < num_layer; i++)
	{
		c_half[i] = (1 / kappa_half[i])*(exp(-alpha_half[i] * source.dt) - 1);
	}
}

void cvl::set_half_coe(src source)
{
	int i;
	cvl_half_coe = (float *)malloc(sizeof(float) * num_layer);
	for (i = 0; i < num_layer; i++)
	{
		cvl_half_coe[i] = exp(-alpha_half[i] * source.dt);
	}
}

void cvl::set_full_coe(src source)
{
	int i;
	cvl_full_coe = (float *)malloc(sizeof(float) * num_layer);
	for (i = 0; i < num_layer; i++)
	{
		cvl_full_coe[i] = exp(-alpha_full[i] * source.dt);
	}
}

void cvl::cmp_cvlh(src s, E Ez, H Hy, int time)
{
	int i, j;
	int pmlbd = num_layer - 1;
	int iopen_side, iclose_side, jopen_side, jclose_side;
	int iopen_ud, iclose_ud, jopen_ud, jclose_ud;
	int offset_e;
	//left & right CPML
	offset_e = num_layer + s.size_x + 1;//51-10=41
	iopen_side = 0; iclose_side = side_sy;
	jopen_side = 0; jclose_side = side_sx;
	for (i = iopen_side; i < iclose_side; i++)
	{
		for (j = jopen_side; j < jclose_side; j++)
		{
			//Hxzl
			Hxzl[i*side_sx + j] = cvl_full_coe[pmlbd - j] * Hxzl[i*side_sx + j]
				+ c_full[pmlbd - j] * (Ez.Ez[(i + 1)*Ez.size_x + j] - Ez.Ez[i*Ez.size_x + j]) / s.dz;

			//Hxzr
			Hxzr[i*side_sx + j] = cvl_full_coe[j] * Hxzr[i*side_sx + j]
				+ c_full[j] * (Ez.Ez[(i + 1)*Ez.size_x + j + offset_e] - Ez.Ez[i*Ez.size_x + j + offset_e]) / s.dz;
		}
	}

	//3, 17
	offset_e = num_layer + s.size_x;
	iopen_side = 0; iclose_side = side_sy + 1;
	jopen_side = 0; jclose_side = side_sx;
	for (i = iopen_side; i < iclose_side; i++)
	{
		for (j = jopen_side; j < jclose_side; j++)
		{
			//Hyzl
			Hyzl[i*side_sx + j] = cvl_half_coe[pmlbd - j] * Hyzl[i*side_sx + j]
				+ c_half[pmlbd - j] * (Ez.Ez[i*Ez.size_x + j + 1] - Ez.Ez[i*Ez.size_x + j]) / s.dz;
			//Hyzr
			Hyzr[i*side_sx + j] = cvl_half_coe[j] * Hyzr[i*side_sx + j]
				+ c_half[j] * (Ez.Ez[i*Ez.size_x + j + 1 + offset_e] - Ez.Ez[i*Ez.size_x + j + offset_e]) / s.dz;
			//			if (time == 6 && ( i == 8))
			//			{
			//				cout << "Hyzr(" << i << ", " << j << "): " << Hyzr[i*side_sx + j] << endl;
			//				cout << "\tDebug\t"
			//					<<  Ez.Ez[i*Ez.size_x + j + 1 + offset_e]<<"\t"
			//					<<Ez.Ez[i*Ez.size_x + j + offset_e]
			//					<< endl;
			//			}
		}
	}

	//bottom & upper CPML, 11, 3
	iopen_ud = 0; iclose_ud = ud_sy;
	jopen_ud = 0; jclose_ud = ud_sx + 1;
	for (i = iopen_ud; i < iclose_ud; i++)
	{
		for (j = jopen_ud; j < jclose_ud; j++)
		{
			//bottom
			offset_e = num_layer;
			//Hxzd
			Hxzd[i*(ud_sx + 1) + j] = cvl_half_coe[pmlbd - i] * Hxzd[i*(ud_sx + 1) + j]
				+ c_half[pmlbd - i] * (Ez.Ez[(i + 1)*Ez.size_x + j + offset_e] - Ez.Ez[i*Ez.size_x + j + offset_e]) / s.dz;

			//upper
			offset_e = num_layer + s.size_y;
			//Hxzu
			Hxzu[i*(ud_sx + 1) + j] = cvl_half_coe[i] * Hxzu[i*(ud_sx + 1) + j]
				+ c_half[i] * (Ez.Ez[(i + 1 + offset_e)*Ez.size_x + j + num_layer] - Ez.Ez[(i + offset_e)*Ez.size_x + j + num_layer]) / s.dz;
		}
	}

	//10, 3
	iopen_ud = 0; iclose_ud = ud_sy;
	jopen_ud = 0; jclose_ud = ud_sx;
	for (i = iopen_ud; i < iclose_ud; i++)
	{
		for (j = jopen_ud; j < jclose_ud; j++)
		{
			//Hyzd
			Hyzd[i*ud_sx + j] = cvl_full_coe[pmlbd - i] * Hyzd[i*ud_sx + j]
				+ c_full[pmlbd - i] * (Ez.Ez[i*Ez.size_x + j + 1 + offset_e] - Ez.Ez[i*Ez.size_x + j + offset_e]) / s.dz;

			//Hyzu
			Hyzu[i*ud_sx + j] = cvl_full_coe[i] * Hyzu[i*ud_sx + j]
				+ c_full[i] * (Ez.Ez[(i + offset_e)*Ez.size_x + j + 1 + num_layer] - Ez.Ez[(i + offset_e)*Ez.size_x + j + num_layer]) / s.dz;
		}
	}
}

void cvl::cmp_cvle(src s, E Ez, H h, int time)
{
	int i, j;
	int pmlbd = num_layer - 1;
	int iopen_side = 0, iclose_side = side_sy;
	int jopen_side = 0, jclose_side = side_sx - 1;
	int iopen_ud = 0, iclose_ud = ud_sy - 1;
	int jopen_ud = 0, jclose_ud = ud_sx;
	int offset_h;

	//left side
	for (i = iopen_side; i <= iclose_side; i++)
	{
		for (j = jopen_side; j <= jclose_side; j++)
		{
			if (i == 0 || j == 0 || i == iclose_side)
			{
				Exyl[i*side_sx + j] = 0.f;
				Eyxl[i*side_sx + j] = 0.f;
				continue;
			}
			//Exyl
			Exyl[i*side_sx + j] = cvl_full_coe[pmlbd - j] * Exyl[i*side_sx + j]
				+ cvl_full_coe[pmlbd - j] * (h.Hy[i*h.size_Hy_x + j] - h.Hy[(i - 1)*h.size_Hy_x + j]) / s.dz;
			//Eyxl
			Eyxl[i*side_sx + j] = cvl_full_coe[pmlbd - j] * Eyxl[i*side_sx + j]
				+ cvl_full_coe[pmlbd - j] * (h.Hx[i*h.size_Hx_x + j] - h.Hx[i*h.size_Hx_x + j - 1]) / s.dz;
			//			if (time == 5 && (i == 8) && j == jclose_side)
			//			{
			//				cout << "Eyxl(" << i << ", " << j << "): " << Eyxl[i*side_sx + j] << endl;
			//				cout << "\tDebug\t"
			//					<< "(" << i << ", " << j << ")" << "\t"
			//					<< h.Hx[i*h.size_Hx_x + j] << "\t"
			//					<< "Hx location: (" << i << ", " << j - 1 << ")\t"
			//					<< h.Hx[i*h.size_Hx_x + j - 1]
			//					<< endl;
			//			}
		}
	}

	//right side
	offset_h = num_layer + s.size_x + 1;
	if (time == 0)
	{
		cout << "iopen_ud: " << iopen_ud << "\ticlose_ud: " << iclose_ud << endl;
		cout << "jopen_ud: " << jopen_ud << "\tjclose_ud: " << jclose_ud << endl;
		cout << "offset_h: " << offset_h << endl;
	}
	for (i = iopen_side; i <= iclose_side; i++)
	{
		for (j = jopen_side; j <= jclose_side; j++)
		{
			if (i == 0 || i == iclose_side || j == jclose_side)
			{
				Exyr[i*side_sx + j] = 0.f;
				Eyxr[i*side_sx + j] = 0.f;
				continue;
			}
			//Exyr
			Exyr[i*side_sx + j] = cvl_full_coe[j] * Exyr[i*side_sx + j]
				+ cvl_full_coe[j] * (h.Hy[i*h.size_Hy_x + j + offset_h] - h.Hy[(i - 1)*h.size_Hy_x + j + offset_h]) / s.dz;
			//Eyxr
			Eyxr[i*side_sx + j] = cvl_full_coe[j] * Eyxr[i*side_sx + j]
				+ cvl_full_coe[j] * (h.Hx[i*h.size_Hx_x + j + offset_h] - h.Hx[i*h.size_Hx_x + j - 1 + offset_h]) / s.dz;
		}
	}

	//upper, (11, 3)
	int width = ud_sx + 1;
	offset_h = num_layer + s.size_y + 1;
	for (i = iopen_ud; i <= iclose_ud; i++)
	{
		for (j = jopen_ud; j <= jclose_ud; j++)
		{
			if (i == iclose_ud)
			{
				Exyu[i*width + j] = 0.f;
				Eyxu[i*width + j] = 0.f;
				continue;
			}
			//Exyu
			Exyu[i*width + j] = cvl_full_coe[i] * Exyu[i*width + j]
				+ cvl_full_coe[i] * (h.Hy[(i + offset_h)*h.size_Hy_x + j + num_layer] - h.Hy[(i - 1 + offset_h)*h.size_Hy_x + j + num_layer]) / s.dz;
			//Eyxu
			Eyxu[i*width + j] = cvl_full_coe[i] * Eyxu[i*width + j]
				+ cvl_full_coe[i] * (h.Hx[(i + offset_h)*h.size_Hx_x + j + num_layer] - h.Hx[(i + offset_h)*h.size_Hx_x + j - 1 + num_layer]) / s.dz;
		}
	}
	//bottom, (11, 3)
	for (i = iopen_ud; i <= iclose_ud; i++)
	{
		for (j = jopen_ud; j <= jclose_ud; j++)
		{
			if (i == 0)
			{
				Exyd[i*width + j] = 0.f;
				Eyxd[i*width + j] = 0.f;
				continue;
			}
			//Exyd
			Exyd[i*width + j] = cvl_full_coe[pmlbd - i] * Exyd[i*width + j]
				+ cvl_full_coe[pmlbd - i] * (h.Hy[i*h.size_Hy_x + j + num_layer] - h.Hy[(i - 1)*h.size_Hy_x + j + num_layer]) / s.dz;
			//Eyxd
			Eyxd[i*width + j] = cvl_full_coe[pmlbd - i] * Eyxd[i*width + j]
				+ cvl_full_coe[pmlbd - i] * (h.Hy[i*h.size_Hy_x + j+ + num_layer] - h.Hy[i*h.size_Hy_x + j + num_layer]) / s.dz;
			if (time == 6 && i == iclose_ud && (j == 6))
			{
				cout << "Eyxd(" << i << ", " << j << "): " << Eyxd[i*width + j] << endl;
				cout << "\tDebug\t"
					<< "width: " << width << "\t"
					<< h.Hy[i*h.size_Hy_x + j+1 + num_layer] << "\t"
					<< "Hy location: (" << i << ", " << j  + num_layer << ")\t"
					<< "Hy number: " << i*h.size_Hy_x + j  + num_layer << "\t"
					<< h.Hy[i*h.size_Hy_x + j + num_layer]
					<< endl;
			}
		}
	}
}

void cvl::checkout()
{
	//for coefficient
	cout << "convolution class checkout:" << endl;
	cout << "\tm = " << m << endl;
	cout << "\teps_r = " << eps_r << endl;
	cout << "\teps0 = " << eps0 << endl;
	cout << "\tsigma_max = " << sigma_max << endl;
	cout << "\tkappa_max = " << kappa_max << endl;
	cout << "\td = " << d << endl;
	cout << endl;

	cout << "\tside grid number x(side_sx):\t" << side_sx << endl;
	cout << "\tside grid number y(side_sy):\t" << side_sy << endl;
	cout << "\tup & down grid number x(ud_sx):\t" << ud_sx << endl;
	cout << "\tup & down grid number y(ud_sy):\t" << ud_sy << endl;
	cout << endl;
	//for variables
	int i;
	//cout << "c_full" << endl;
	//for (i = 0; i < num_layer; i++)
	//{
	//	cout << c_full[i] << endl;
	//}
}

void cvl::save2file()
{
	save2file_e();
	save2file_h();
}

void cvl::save2file_h()
{
	fstream myfile;
	int i, j;

	myfile.open("Hxzl.txt", ios::app);
	for (i = side_sy - 1; i >= 0; i--)
	{
		for (j = 0; j < side_sx; j++)
		{
			myfile << Hxzl[i*side_sx + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();

	myfile.open("Hyzl.txt", ios::app);
	for (i = side_sy; i >= 0; i--)
	{
		for (j = 0; j < side_sx; j++)
		{
			myfile << Hyzl[i*side_sx + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();

	myfile.open("Hyzr.txt", ios::app);
	for (i = side_sy; i >= 0; i--)
	{
		for (j = 0; j < side_sx; j++)
		{
			myfile << Hyzr[i*side_sx + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();

	myfile.open("Hyzu.txt", ios::app);
	for (i = ud_sy - 1; i >= 0; i--)
	{
		for (j = 0; j < ud_sx; j++)
		{
			myfile << Hyzu[i*ud_sx + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();

	myfile.open("Hxzu.txt", ios::app);
	for (i = ud_sy - 1; i >= 0; i--)
	{
		for (j = 0; j < ud_sx + 1; j++)
		{
			myfile << Hxzu[i*(ud_sx + 1) + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();

}
void cvl::save2file_e()
{
	fstream myfile;
	int i, j;

	myfile.open("Exyl.txt", ios::app);
	for (i = side_sy; i >= 0; i--)
	{
		for (j = 0; j < side_sx; j++)
		{
			myfile << Exyl[i*side_sx + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();

	myfile.open("Exyr.txt", ios::app);
	for (i = side_sy; i >= 0; i--)
	{
		for (j = 0; j < side_sx; j++)
		{
			myfile << Exyr[i*side_sx + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();

	myfile.open("Eyxl.txt", ios::app);
	for (i = side_sy; i >= 0; i--)
	{
		for (j = 0; j < side_sx; j++)
		{
			myfile << Eyxl[i*side_sx + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();

	myfile.open("Eyxr.txt", ios::app);
	for (i = side_sy; i >= 0; i--)
	{
		for (j = 0; j < side_sx; j++)
		{
			myfile << Eyxr[i*side_sx + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();

}

void cvl::save2file_coe()
{
	//fstream myfile;
	//int i;

	//myfile.open("coe.txt", ios::app);
	//for (i = 0; i < num_layer; i++)
	//{
	//	myfile << i << "\t";
	//	myfile << "distance= " << distance_half[i] << "\t";
	//	myfile << "kappa= " << kappa_half[i] << "\t";
	//	myfile << "sigma= " << sigma_half[i] << "\t";
	//	myfile << "alpha= " << alpha_half[i] << "\t";
	//	myfile << endl;
	//}
	//myfile << endl;
	//for (i = 0; i < num_layer; i++)
	//{
	//	myfile << i << "\t";
	//	myfile << "distance= " << distance_full[i] << "\t";
	//	myfile << "kappa= " << kappa_full[i] << "\t";
	//	myfile << "sigma= " << sigma_full[i] << "\t";
	//	myfile << "alpha= " << alpha_full[i] << "\t";
	//	myfile << endl;
	//}
	//myfile.close();
}