#include "E.h"
#define MUR

E::E(src s)
{
	Ex_nbdl = 0.f;
	Ex_nbdr = 0.f;
	Ex_bdl = 0.f;
	Ex_bdr = 0.f;

	size_Ex = s.size_space + 1;
	Ex = (float *)malloc(size_Ex * sizeof(float));
	memset(Ex, 0, size_Ex * sizeof(float));

	coe_E = s.dt / (epsilon * s.dz);
	coe_E_cvl = s.dt / epsilon;
	coe_mur = (C*s.dt - s.dz) / (C*s.dt + s.dz);

	fstream myfile;
	myfile.open("Ex.txt", ios::out);
	myfile.close();
}

void E::cmp(H Hy)
{
	int i;
	for (i = 1; i < size_Ex-1; i++){
		Ex[i] = Ex[i] - (coe_E)*(Hy.Hy[i] - Hy.Hy[i - 1]);
	}
}

void E::cmp(H Hy, cvl cvln)
{
	int i;
	for (i = 1; i < size_Ex-1; i++){
		Ex[i] = Ex[i] - (coe_E)*(Hy.Hy[i] - Hy.Hy[i - 1]);
	}
	Ex[0] = Ex[0] - coe_E*(Hy.Hy[0] - cvln.Hyl[cvln.num_layer - 1]);
	Ex[size_Ex - 1] = Ex[size_Ex - 1] - coe_E*(cvln.Hyr[0] - Hy.Hy[Hy.size_Hy - 1]);
}

void E::checkout()
{
	cout << "Ex: " << endl;
	for (int i=0; i < size_Ex; i++){
		cout << Ex[i] << "\t"; 
	}
	cout << endl;
}

void E::boundary()
{
#ifndef MUR
	Ex[0] = 0.f;
	Ex[size_Ex - 1] = 0.f;
#else
	int bdr, nbdr, bdl, nbdl;
	bdr = size_Ex - 1;
	nbdr = bdr - 1;
	bdl = 0;
	nbdl = 1;

	Ex[bdr] = Ex_nbdr + coe_mur*(Ex[nbdr] - Ex_bdr);
	Ex[bdl] = Ex_nbdl + coe_mur*(Ex[nbdl] - Ex_bdl);

	Ex_bdr = Ex[bdr];
	Ex_nbdr = Ex[nbdr];
	Ex_bdl = Ex[bdl];
	Ex_nbdl = Ex[nbdl];
#endif
}


void E::save2file()
{
	fstream myfile;
	myfile.open("Ex.txt", ios::app);

	for (int i = 0; i < size_Ex; i++)
	{
		myfile << Ex[i] << "\t";
	}
	myfile << endl;
	myfile.close();
}