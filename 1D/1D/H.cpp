#include "H.h"

H::H(src source)
{
	int i;
	size_Hy = source.size_space;

	Hy = (float *)malloc(size_Hy * sizeof(float));
	memset(Hy, 0, size_Hy * sizeof(float));

	coe_H = source.dt / (mu * source.dz);

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