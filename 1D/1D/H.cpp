#include "H.h"

H::H(src s)
{
	int i;
	size_Hy = s.size_space;

	Hy = (float *)malloc(size_Hy * sizeof(float));
	memset(Hy, 0, size_Hy * sizeof(float));

	coe_H = s.dt / (mu * s.dz);
	coe_H_cvl = s.dt / mu;

	fstream myfile;
	myfile.open("Hy.txt", ios::out);
	myfile.close();
}

void H::cmp(E Ex)
{
	int i;
	for (i = 0; i < size_Hy; i++){
		Hy[i] = Hy[i] - (coe_H)*(Ex.Ex[i + 1] - Ex.Ex[i]);
	}
}

void H::checkout()
{
	cout << "Hy: " << endl;
	for (int i=0; i < size_Hy; i++){
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