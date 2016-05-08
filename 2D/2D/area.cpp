#include "area.h"

area::area(int L, string name)
{
	length = L;
	grid_num = length;
	filename = name;
	p.insert(p.end(), grid_num, 0);
	itb = p.begin();
	ite = p.end();

	fstream my;
	my.open(filename, ios::out);
	my.close();
}


area::~area()
{
}

void area::checkout()
{
	cout << filename << endl;
	cout << "\tlength: " << length << endl;
	if (!p.empty())
	{
		cout << "\tp is not empty" << endl;
		cout << "\tthe first element of p: " << p.front() << endl;
		cout << "\tthe last element of p: " << p.back() << endl;
	}
}

//void area::real_pos(int i, int j)
//{
//}

void area::save2file(string name)
{
	int j;
	fstream myfile;
	myfile.open(name, ios::app);
	for (j = 0; j < length; j++){
		myfile << p[j] << "\t";
	}
	myfile << endl;
	myfile.close();
}