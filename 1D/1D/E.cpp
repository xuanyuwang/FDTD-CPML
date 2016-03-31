#include "E.h"

E::E(src source)
{
	int i;

	size_Ex = source.size_space + 1;
	Ex = (float *)malloc(size_Ex * sizeof(float));
	memset(Ex, 0, size_Ex * sizeof(float));

	coe_E = source.dt / (epsilon * source.dz);

	fstream myfile;
	myfile.open("Ex.txt", ios::out);
	myfile.close();
}

void E::cmp(H Hy)
{
	int i;
	for (i = 1; i < size_Ex; i++){
		Ex[i] = Ex[i] - (coe_E)*(Hy.Hy[i] - Hy.Hy[i - 1]);
	}
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
	Ex[0] = 0.f;
	Ex[size_Ex - 1] = 0.f;
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