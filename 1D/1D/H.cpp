#include "H.h"
//#define debug

H::H(src s, cvl c)
{
	int i;
	size_Hy = s.size_space + 2 * c.num_layer;

	Hy = (float *)malloc(size_Hy * sizeof(float));
	memset(Hy, 0, size_Hy * sizeof(float));

	coe_H = s.dt / (mu * s.dz);
	coe_H_cvl = s.dt / mu;

	fstream myfile;
	myfile.open("Hy.txt", ios::out);
	for (i = 0; i < size_Hy; i++)
	{
		myfile << i << "\t";
	}
	myfile << endl;
	myfile.close();
}

void H::cmp(E Ex, cvl c, src s, int time)
{
	int i;
	int pmlbd = c.num_layer - 1;
	int sec1l = 0, sec1r = pmlbd;//[0, 9]
	int sec2l = pmlbd + 1, sec2r = pmlbd + s.size_space;//[10, 39]
	int sec3l = size_Hy - c.num_layer, sec3r = size_Hy - 1;//[40, 49]
	//left CPML. [0, pmlbd]
	for (i = sec1l; i <= sec1r; i++)
	{
		Hy[i] = Hy[i] - (s.dt / (mu*c.kappa_H[pmlbd - i] * s.dz))*(Ex.Ex[i + 1] - Ex.Ex[i])
			- (s.dt / mu)*(c.Hzxl[i]);
	}
	//FDTD area
	for (i = sec2l; i <= sec2r; i++)
	{
		Hy[i] = Hy[i] - coe_H*(Ex.Ex[i + 1] - Ex.Ex[i]);
	}
	//right CPML
	for (i = sec3l; i <= sec3r; i++)
	{
		Hy[i] = Hy[i] - (s.dt / (mu*c.kappa_H[i - sec3l] * s.dz))*(Ex.Ex[i + 1] - Ex.Ex[i])
			- (s.dt / mu)*(c.Hzxr[i - sec3l]);
	}
#ifdef debug
	if (time > 12 && time < 16)
	{
		cout << time << ": ";
		cout << Hy[sec3l - 2] << "\t" << Hy[sec3l - 1] << endl;
	}
	if (time == 15)
	{
		cout << sec1l << "\t" << sec1r << endl;
		cout << sec2l << "\t" << sec2r << endl;
		cout << sec3l << "\t" << sec3r << endl;
		cout << (Hy[sec2r] - Hy[sec2r - 1]) << endl;
	}
#endif
}

void H::checkout()
{
	cout << "Hy: " << endl;
	for (int i = 0; i < size_Hy; i++){
		cout << Hy[i] << "\t";
	}
	cout << endl;
}

void H::save2file()
{
	fstream myfile;
	myfile.open("Hy.txt", ios::app);

	for (int i = 0; i < size_Hy; i++)
	{
		myfile << Hy[i] << "\t";
	}
	myfile << endl;
	myfile.close();
}