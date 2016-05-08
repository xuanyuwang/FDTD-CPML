#include "area.h"

area::area(int w, int h, string name)
{
	width = w;
	height = h;
	grid_num = width*height;
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
	cout << "\twidth: " << width << endl;
	cout << "\theight: " << height << endl;
	if (!p.empty())
	{
		cout << "\tp is not empty" << endl;
		cout << "\tthe first element of p: " << p.front() << endl;
		cout << "\tthe last element of p: " << p.back() << endl;
	}
}

void area::save2file(string name)
{
	int i, j;
	fstream myfile;
	myfile.open(name, ios::app);
	for (i = height - 1; i >= 0; i--)
	{
		for (j = 0; j < width; j++)
		{
			myfile << p.at(i*width + j) << "\t";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}
