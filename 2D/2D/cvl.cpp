#include "cvl.h"
#define debug

/*Class COE*/
COE::COE(int number_layer, src s)
{
	num_layer = number_layer;
	d = number_layer*s.dz;
	sigma_max = (m + 1) / (sqrt(eps_r) * 150 * PI * s.dz);
	side_sx = num_layer;
	side_sy = s.size_y + 2 * num_layer;
	ud_sx = s.size_x;//10
	ud_sy = num_layer;//3

	set_distance(s);
	set_kappa(s);
	set_sigma(s);
	set_alpha(s);
	set_c(s);
	set_full_coe(s);
	set_half_coe(s);

	fstream myfile;
	myfile.open("coe.txt", ios::out);
	myfile.close();
}

//the distance of every field points of H and E in CPML layer.
//the "i_min" is the nearest distance, the "I_max" is the furthest distance.
void COE::set_distance(src source)
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

void COE::set_sigma(src source)
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

void COE::set_kappa(src source)
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

void COE::set_alpha(src source)
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

void COE::set_c(src source)
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

void COE::set_half_coe(src source)
{
	int i;
	cvl_half_coe = (float *)malloc(sizeof(float) * num_layer);
	for (i = 0; i < num_layer; i++)
	{
		cvl_half_coe[i] = exp(-alpha_half[i] * source.dt);
	}
}

void COE::set_full_coe(src source)
{
	int i;
	cvl_full_coe = (float *)malloc(sizeof(float) * num_layer);
	for (i = 0; i < num_layer; i++)
	{
		cvl_full_coe[i] = exp(-alpha_full[i] * source.dt);
	}
}

void COE::checkout()
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

COE::~COE()
{
	//free(distance_half);
	//free(distance_full);
	//free(sigma_half);
	//free(sigma_full);
	//free(alpha_full);
	//free(alpha_half);
	//free(kappa_half);
	//free(kappa_full);
	//free(cvl_half_coe);
	//free(cvl_full_coe);
	//free(c_full);
	//free(c_half);
}

/**************************** Hxz*******************************/
/*Class HXZL*/
HXZL::HXZL(COE c, src s)
{
	width = c.num_layer;
	height = s.size_y + 2 * c.num_layer;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Hxzl.txt", ios::out);
	myfile.close();
}

void HXZL::cmp(COE c, src s, E e, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			p[i*width + j] = c.cvl_full_coe[pmlbd - j] * p[i*width + j]
				+ c.c_full[pmlbd - j]
				* (e.Ez[(i + 1)*e.size_x + j] - e.Ez[i*e.size_x + j]) / s.dz;
		}
	}
}

void HXZL::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Hxzl.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*Class HXZR*/
HXZR::HXZR(COE c, src s)
{
	width = c.num_layer;
	height = s.size_y + 2 * c.num_layer;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Hxzr.txt", ios::out);
	myfile.close();
}

void HXZR::cmp(COE c, src s, E e, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;
	int offset_e = c.num_layer + s.size_x + 1;//

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			p[i*width + j] = c.cvl_full_coe[j] * p[i*width + j]
				+ c.c_full[j] * (e.Ez[(i + 1)*e.size_x + j + offset_e] - e.Ez[i*e.size_x + j + offset_e]) / s.dz;
		}
	}
}

void HXZR::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Hxzr.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*Class HXZU*/
HXZU::HXZU(COE c, src s)
{
	width = s.size_x + 1;
	height = c.num_layer;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Hxzu.txt", ios::out);
	myfile.close();
}

void HXZU::cmp(COE c, src s, E e, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;
	int offset_e = c.num_layer + s.size_y;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			p[i*width + j] = c.cvl_half_coe[i] * p[i*width + j]
				+ c.c_half[i] * (e.Ez[(i + 1 + offset_e)*e.size_x + j + c.num_layer]
				- e.Ez[(i + offset_e)*e.size_x + j + c.num_layer]) / s.dz;
		}
	}
}

void HXZU::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Hxzu.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*Class HXZD*/
HXZD::HXZD(COE c, src s)
{
	width = s.size_x + 1;
	height = c.num_layer;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Hxzd.txt", ios::out);
	myfile.close();
}

void HXZD::cmp(COE c, src s, E e, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;
	int offset_e = c.num_layer;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			p[i*width + j] = c.cvl_half_coe[pmlbd - i] * p[i*width + j]
				+ c.c_half[pmlbd - i] * (e.Ez[(i + 1)*e.size_x + j + offset_e]
				- e.Ez[i*e.size_x + j + offset_e]) / s.dz;
		}
	}
}

void HXZD::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Hxzd.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*************************Hyz*******************************/
/*Class HYZL*/
HYZL::HYZL(COE c, src s)
{
	width = c.num_layer;
	height = s.size_y + 2 * c.num_layer + 1;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Hyzl.txt", ios::out);
	myfile.close();
}

void HYZL::cmp(COE c, src s, E e, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			p[i*width + j] = c.cvl_half_coe[pmlbd - j] * p[i*width + j]
				+ c.c_half[pmlbd - j] * (e.Ez[i*e.size_x + j + 1]
				- e.Ez[i*e.size_x + j]) / s.dz;
		}
	}
}

void HYZL::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Hyzl.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*Class HYZR*/
HYZR::HYZR(COE c, src s)
{
	width = c.num_layer;
	height = s.size_y + 2 * c.num_layer + 1;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Hyzr.txt", ios::out);
	myfile.close();
}

void HYZR::cmp(COE c, src s, E e, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;
	int offset_e = c.num_layer + s.size_x;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			p[i*width + j] = c.cvl_half_coe[j] * p[i*width + j]
				+ c.c_half[j] * (e.Ez[i*e.size_x + j + 1 + offset_e]
				- e.Ez[i*e.size_x + j + offset_e]) / s.dz;
		}
	}
}

void HYZR::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Hyzr.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*Class HYZU*/
HYZU::HYZU(COE c, src s)
{
	width = s.size_x;
	height = c.num_layer;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Hyzu.txt", ios::out);
	myfile.close();
}

void HYZU::cmp(COE c, src s, E e, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;
	int offset_e = c.num_layer + s.size_y + 1;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			p[i*width + j] = c.cvl_full_coe[i] * p[i*width + j]
				+ c.c_full[i] * (e.Ez[(i + offset_e)*e.size_x + j + 1 + c.num_layer]
				- e.Ez[(i + offset_e)*e.size_x + j + c.num_layer]) / s.dz;
		}
	}
}

void HYZU::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Hyzu.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*Class HYZD*/
HYZD::HYZD(COE c, src s)
{
	width = s.size_x + 1;
	height = c.num_layer;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Hyzd.txt", ios::out);
	myfile.close();
}

void HYZD::cmp(COE c, src s, E e, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;
	int offset_e = c.num_layer;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			p[i*width + j] = c.cvl_full_coe[pmlbd - i] * p[i*width + j]
				+ c.c_full[pmlbd - i] * (e.Ez[i*e.size_x + j + 1 + offset_e]
				- e.Ez[i*e.size_x + j + offset_e]) / s.dz;
		}
	}
}

void HYZD::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Hyzd.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*************************Exy*******************************/
/*Class EXYL*/
EXYL::EXYL(COE c, src s)
{
	width = c.num_layer;
	height = s.size_y + 2 * c.num_layer + 1;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Exyl.txt", ios::out);
	myfile.close();
}

void EXYL::cmp(COE c, src s, HY hy, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			if (i == 0 || j == 0 || i == height - 1)
			{
				p[i*width + j] = 0.f;
				continue;
			}
			p[i*width + j] = c.cvl_full_coe[pmlbd - j] * p[i*width + j]
				+ c.c_full[pmlbd - j] * (hy.p[i*hy.width + j]
				- hy.p[(i - 1)*hy.width + j]) / s.dz;
		}
	}
}

void EXYL::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Exyl.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*Class ExyR*/
EXYR::EXYR(COE c, src s)
{
	width = c.num_layer;
	height = s.size_y + 2 * c.num_layer + 1;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Exyr.txt", ios::out);
	myfile.close();
}

void EXYR::cmp(COE c, src s,HY hy, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;
	int offset_h = c.num_layer + s.size_x + 1;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			if (i == 0 || i == height - 1 || j == width - 1)
			{
				p[i*width + j] = 0.f;
				continue;
			}
			p[i*width + j] = c.cvl_full_coe[pmlbd - j] * p[i*width + j]
				+ c.c_full[pmlbd - j] * (hy.p[i*hy.width + j + offset_h]
				- hy.p[(i - 1)*hy.width + j + offset_h]) / s.dz;
		}
	}
}

void EXYR::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Exyr.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*Class EXYU*/
EXYU::EXYU(COE c, src s)
{
	width = s.size_x + 1;
	height = c.num_layer;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Exyu.txt", ios::out);
	myfile.close();
}

void EXYU::cmp(COE c, src s,HY hy, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;
	int offset_h = c.num_layer + s.size_y + 1;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			if (i == height - 1)
			{
				p[i*width + j] = 0.f;
				continue;
			}
			p[i*width + j] = c.cvl_full_coe[pmlbd - j] * p[i*width + j]
				+ c.c_full[pmlbd - j] * (hy.p[(i + offset_h)*hy.width + j + c.num_layer]
				- hy.p[(i - 1 + offset_h)*hy.width + j + c.num_layer]) / s.dz;
		}
	}
}

void EXYU::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Hyzu.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*Class EXYD*/
EXYD::EXYD(COE c, src s)
{
	width = s.size_x + 1;
	height = c.num_layer;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Exyd.txt", ios::out);
	myfile.close();
}

void EXYD::cmp(COE c, src s,HY hy, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			if (i == 0)
			{
				p[i*width + j] = 0.f;
				continue;
			}
			p[i*width + j] = c.cvl_full_coe[pmlbd - j] * p[i*width + j]
				+ c.c_full[pmlbd - j] * (hy.p[i*hy.width + j + c.num_layer]
				- hy.p[(i - 1)*hy.width + j + c.num_layer]) / s.dz;
		}
	}
}

void EXYD::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Exyd.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*************************Eyx*******************************/
/*Class EYXL*/
EYXL::EYXL(COE c, src s)
{
	width = c.num_layer;
	height = s.size_y + 2 * c.num_layer + 1;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Eyxl.txt", ios::out);
	myfile.close();
}

void  EYXL::cmp(COE c, src s,HX hx, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			if (i == 0 || j == 0 || i == height - 1)
			{
				p[i*width + j] = 0.f;
				continue;
			}
			p[i*width + j] = c.cvl_full_coe[pmlbd - j] * p[i*width + j]
				+ c.c_full[pmlbd - j] * (hx.p[i*hx.width + j] - hx.p[i*hx.width + j - 1]) / s.dz;
		}
	}
}

void EYXL::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Eyxl.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*Class EyxR*/
EYXR::EYXR(COE c, src s)
{
	width = c.num_layer;
	height = s.size_y + 2 * c.num_layer + 1;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Eyxr.txt", ios::out);
	myfile.close();
}

void EYXR::cmp(COE c, src s,HX hx, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;
	int offset_h = c.num_layer + s.size_x + 1;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			if (i == 0 || i == height - 1 || j == width - 1)
			{
				p[i*width + j] = 0.f;
				continue;
			}
			p[i*width + j] = c.cvl_full_coe[pmlbd - j] * p[i*width + j]
				+ c.c_full[pmlbd - j] * (hx.p[i*hx.width + j + offset_h]
				- hx.p[i*hx.width + j - 1 + offset_h]) / s.dz;
		}
	}
}

void EYXR::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Eyxr.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*Class EYXU*/
EYXU::EYXU(COE c, src s)
{
	width = s.size_x + 1;
	height = c.num_layer;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Eyxu.txt", ios::out);
	myfile.close();
}

void EYXU::cmp(COE c, src s,HX hx, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;
	int offset_h = c.num_layer + s.size_y + 1;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			if (i == height - 1)
			{
				p[i*width + j] = 0.f;
				continue;
			}
			p[i*width + j] = c.cvl_full_coe[pmlbd - j] * p[i*width + j]
				+ c.c_full[pmlbd - j] * (hx.p[(i + offset_h)*hx.width + j + c.num_layer]
				- hx.p[(i + offset_h)*hx.width + j - 1 + c.num_layer]) / s.dz;
		}
	}
}

void EYXU::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Hyzu.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

/*Class EYXD*/
EYXD::EYXD(COE c, src s)
{
	width = s.size_x + 1;
	height = c.num_layer;

	TOOLS::alloc(p, width, height);

	fstream myfile;
	myfile.open("Eyxd.txt", ios::out);
	myfile.close();
}

void EYXD::cmp(COE c, src s,HX hx, int time)
{
	int i, j;
	int pmlbd = c.num_layer - 1;

	for (i = 0; i < height; i++)
	{
		for (j = 0; j < width; j++)
		{
			if (i == 0)
			{
				p[i*width + j] = 0.f;
				continue;
			}
			p[i*width + j] = c.cvl_full_coe[pmlbd - j] * p[i*width + j]
				+ c.c_full[pmlbd - j] * (hx.p[i*hx.width + j + 1 + c.num_layer]
				- hx.p[i*hx.width + j + c.num_layer]) / s.dz;
		}
	}
}

void EYXD::save2file()
{
	int i, j;
	fstream myfile;

	myfile.open("Eyxd.txt", ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p[i*width + j] << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}
