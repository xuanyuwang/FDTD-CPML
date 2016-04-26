#include "cvl.h"
#define debug

cvl::cvl(int number, src s)
{
	num_layer = number;
	d = num_layer * s.dz;
	sigma_max = (m + 1) / (sqrt(eps_r) * 150 * PI * s.dz);
	side_sx = num_layer;
	side_sy = s.size_y + 2 * num_layer;
	ud_sx = s.size_x;
	ud_sy = num_layer;

	alloc(Hzxl, side_sx, side_sy);
	alloc(Hzxr, side_sx, side_sy);
	alloc(Ezyl, side_sx, side_sy + 1);
	alloc(Ezyr, side_sx, side_sy + 1);
	alloc(Hzxu, ud_sx, num_layer);
	alloc(Hzxd, ud_sx, num_layer);
	alloc(Ezyu, ud_sx + 1, num_layer);
	alloc(Ezyd, ud_sx + 1, num_layer);

	set_distance(s);
	set_kappa(s);
	set_sigma(s);
	set_alpha(s);
	set_c(s);
	set_E_coe(s);
	set_H_coe(s);

	int i;
	fstream myfile;
	myfile.open("Hzx.txt", ios::out);
	myfile.close();

	myfile.open("Ezy.txt", ios::out);
	myfile.close();

	myfile.open("coe.txt", ios::out);
	myfile.close();
}


cvl::~cvl()
{
	free(Hzxl);
	free(Hzxr);
	free(Hzxu);
	free(Hzxd);
	free(Ezyl);
	free(Ezyr);
	free(Ezyu);
	free(Ezyd);
	free(distance_E);
	free(distance_H);
	free(sigma_E);
	free(sigma_H);
	free(alpha_E);
	free(alpha_H);
	//free(kappa_E);
	free(kappa_H);
	free(c_E);
	free(c_H);
	free(cvl_H_coe);
	free(cvl_E_coe);
}

void cvl::alloc(float* &p, int x, int y)
{
	p = (float*)malloc(x * y * sizeof(float));
	memset(p, 0, x*y*sizeof(float));
}

//the distance of every field points of H and E in CPML layer.
//the "i_min" is the nearest distance, the "I_max" is the furthest distance.
void cvl::set_distance(src source)
{
	int i;
	distance_E = (float *)malloc(num_layer*sizeof(float));
	distance_H = (float *)malloc(num_layer*sizeof(float));
	//for distance_E
	for (i = 0; i < num_layer; i++){
		distance_E[i] = (i + 1) * source.dz;
	}
	//for distance_H
	for (i = 0; i < num_layer; i++){
		distance_H[i] = (i + 0.5) * source.dz;
	}
}

void cvl::set_sigma(src source)
{
	int i;
	sigma_E = (float *)malloc(num_layer * sizeof(float));
	sigma_H = (float *)malloc(num_layer * sizeof(float));
	//sigma E
	for (i = 0; i < num_layer; i++)
	{
		sigma_E[i] = sigma_max * powf((distance_E[i] / d), m);
	}

	//sigma_H
	for (i = 0; i < num_layer; i++)
	{
		sigma_H[i] = sigma_max * powf((distance_H[i] / d), m);
	}
}

void cvl::set_kappa(src source)
{
	kappa_E = (float *)malloc(num_layer * sizeof(float));
	kappa_H = (float *)malloc(num_layer * sizeof(float));
	int i;
	//for kappa_E
	for (i = 0; i <= num_layer; i++)
	{
		kappa_E[i] = 1 + (kappa_max - 1) * powf((distance_E[i] / d), m);
	}
	//for kappa_H
	for (i = 0; i < num_layer; i++)
	{
		kappa_H[i] = 1 + (kappa_max - 1) * powf((distance_H[i] / d), m);
	}
}

void cvl::set_alpha(src source)
{
	alpha_E = (float *)malloc(num_layer * sizeof(float));
	alpha_H = (float *)malloc(num_layer * sizeof(float));

	int i;
	//for alpha_E
	for (i = 0; i < num_layer; i++)
	{
		alpha_E[i] = sigma_E[i] / (eps0 * kappa_E[i]);
	}
	//for alpha_H
	for (i = 0; i < num_layer; i++)
	{
		alpha_H[i] = sigma_H[i] / (eps0 * kappa_H[i]);
	}
}

void cvl::set_c(src source)
{
	c_E = (float *)malloc(num_layer * sizeof(float));
	c_H = (float *)malloc(num_layer * sizeof(float));

	int i;
	//for E
	for (i = 0; i < num_layer; i++)
	{
		c_E[i] = (1 / kappa_E[i])*(exp(-alpha_E[i] * source.dt) - 1);
	}
	//for H
	for (i = 0; i < num_layer; i++)
	{
		c_H[i] = (1 / kappa_H[i])*(exp(-alpha_H[i] * source.dt) - 1);
	}
}

void cvl::set_H_coe(src source)
{
	int i;
	cvl_H_coe = (float *)malloc(sizeof(float) * num_layer);
	for (i = 0; i < num_layer; i++)
	{
		cvl_H_coe[i] = exp(-alpha_H[i] * source.dt);
	}
}

void cvl::set_E_coe(src source)
{
	int i;
	cvl_E_coe = (float *)malloc(sizeof(float) * num_layer);
	for (i = 0; i < num_layer; i++)
	{
		cvl_E_coe[i] = exp(-alpha_E[i] * source.dt);
	}
}

void cvl::cmp_cvlh(src s, E Ex, H Hy, int time)
{
	//int i;
	//int pmlbd = num_layer - 1;
	//int offset_e = Ex.size_Ex - num_layer;//51-10=41
	////Hzxl
	//for (i = 0; i < num_layer; i++)
	//{
	//	Hzxl[i] = cvl_H_coe[pmlbd - i] * Hzxl[i] + c_H[pmlbd - i] * (Ex.Ex[i + 1] - Ex.Ex[i]) / s.dz;
	//}
	////Hzxr
	//for (i = 0; i < num_layer; i++)
	//{
	//	Hzxr[i] = cvl_H_coe[i] * Hzxr[i] + c_H[i] * (Ex.Ex[i + offset_e] - Ex.Ex[i - 1 + offset_e]) / s.dz;
	//}
}

void cvl::cmp_cvle(src s, E ex, H hy, int time)
{
	//int i;
	//int pmlbd = num_layer - 1;
	//int offset_h = hy.size_Hy - num_layer;//40
	////Ezyl
	//Ezyl[0] = Ezyl[1];
	//for (i = 1; i < num_layer; i++)
	//{
	//	Ezyl[i] = cvl_E_coe[pmlbd - i] * Ezyl[i] + c_E[pmlbd - i] * (hy.Hy[i] - hy.Hy[i - 1]) / s.dz;
	//}
	////Ezyr
	//Ezyr[num_layer - 1] = Ezyr[num_layer - 2];
	//for (i = 0; i < num_layer - 1; i++)
	//{
	//	Ezyr[i] = cvl_E_coe[i] * Ezyr[i] + c_E[i] * (hy.Hy[i + 1 + offset_h] - hy.Hy[i + offset_h]) / s.dz;
	//	if (time == 15)
	//	{
	//		cout << i << "\t" << i + offset_h << endl;
	//	}
	//}
}

void cvl::checkout()
{
	//for coefficient
	cout << "m = " << m << endl;
	cout << "eps_r = " << eps_r << endl;
	cout << "eps0 = " << eps0 << endl;
	cout << "sigma_max = " << sigma_max << endl;
	cout << "kappa_max = " << kappa_max << endl;
	cout << "d = " << d << endl;

	//for variables
}

void cvl::save2file()
{
	fstream myfile;
	int i;

	myfile.open("Hzx.txt", ios::app);
	for (i = 0; i < num_layer; i++)
	{
		myfile << Hzxl[i] << "\t";
	}
	for (i = 0; i < num_layer; i++)
	{
		myfile << Hzxr[i] << "\t";
	}
	myfile << endl;
	myfile.close();

	myfile.open("Ezy.txt", ios::app);
	for (i = 0; i < num_layer; i++)
	{
		myfile << Ezyl[i] << "\t";
	}
	for (i = 0; i < num_layer; i++)
	{
		myfile << Ezyr[i] << "\t";
	}
	myfile << endl;
	myfile.close();
}

void cvl::save2file_coe()
{
	fstream myfile;
	int i;

	myfile.open("coe.txt", ios::app);
	for (i = 0; i < num_layer; i++)
	{
		myfile << i << "\t";
		myfile << "distance= " << distance_H[i] << "\t";
		myfile << "kappa= " << kappa_H[i] << "\t";
		myfile << "sigma= " << sigma_H[i] << "\t";
		myfile << "alpha= " << alpha_H[i] << "\t";
		myfile << endl;
	}
	myfile.close();
}