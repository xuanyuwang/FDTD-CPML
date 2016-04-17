#include "E.h"
#define MUR
//#define debug

E::E(src s, cvl c)
{
	size_Ex = (s.size_space + 1) + 2 * c.num_layer;
	Ex = (float *)malloc(size_Ex * sizeof(float));
	memset(Ex, 0, size_Ex * sizeof(float));

	coe_E = s.dt / (epsilon * s.dz);
	coe_E_cvl = s.dt / epsilon;
	//coe_mur = (C*s.dt - s.dz) / (C*s.dt + s.dz);

	fstream myfile;
	myfile.open("Ex.txt", ios::out);
	//for (int i = 0; i < size_Ex; i++)
	//{
		//myfile << i << "\t";
	//}
	myfile << endl;
	myfile.close();
}

void E::cmp(H Hy, cvl c, src s, int time)
{
	int i;
	//[0,pmlbd]([0,9]), [pmlbd,pmlbd+s.size_spzce]([10, 40]), [size_Ex - 10, size_Ex - 1]([41,50])
	int pmlbd = c.num_layer - 1;
	int sec1l = 0, sec1r = pmlbd;
	int sec2l = pmlbd+1, sec2r = pmlbd + s.size_space + 1;
	int sec3l = size_Ex - c.num_layer, sec3r = size_Ex - 1;
	//for left CPML
	for (i = sec1l; i <= sec1r; i++)
	{
		Ex[i] += (s.dt / (epsilon * c.kappa_E[pmlbd - i] * s.dz))*(-(Hy.Hy[i] - Hy.Hy[i - 1]))
			+ (s.dt / epsilon)*(-c.Ezyl[i]);
	}
	//for FDTD area
	for (i = sec2l; i <= sec2r; i++)
	{
		Ex[i] = Ex[i] - coe_E*(Hy.Hy[i] - Hy.Hy[i - 1]);
	}
	//for right CPML
	for (i = sec3l; i <= sec3r; i++)
	{
		Ex[i] += (s.dt / (epsilon*c.kappa_E[i - sec3l] * s.dz))*(-(Hy.Hy[i] - Hy.Hy[i - 1]))
			+ (s.dt / epsilon)*(-c.Ezyr[i - sec3l]);

	}
#ifdef debug
	if (time > 12 && time < 16)
	{
		cout << time << ": ";
		cout << Ex[sec3l - 2] << "\t" << Ex[sec3l - 1] << endl;
	}
	if (time == 15)
	{
		cout << sec1l << "\t" << sec1r << endl;
		cout << sec2l << "\t" << sec2r << endl;
		cout << sec3l << "\t" << sec3r << endl;
		cout << (Hy.Hy[sec2r] - Hy.Hy[sec2r - 1]) << endl;
	}
#endif
	//for boundary
	Ex[0] = 0.f;
	Ex[size_Ex - 1] = 0.f;
}

void E::checkout()
{
	cout << "Ex: " << endl;
	for (int i = 0; i < size_Ex; i++){
		cout << Ex[i] << "\t";
	}
	cout << endl;
}

void E::save2file()
{
	fstream myfile;
	myfile.open("Ex.txt", ios::app);

	for (int i = 10; i <= 40; i++)
	{
		myfile << Ex[i] << "\t";
	}
	myfile << endl;
	myfile.close();
}