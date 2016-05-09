#include "area.h"

/************************************************************************/
/* Name: area::area
 * Arguments: int w, int h, string name.
 * Function: Initialize every member of this class.
 * Description: width = w, height = h, filename = name.
 *	the size of p is is grid_num and all elements are initialize as 0.
 *	Create a file whose name if filename, or clean the content of such a file.
 *	*/
/************************************************************************/
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

/************************************************************************/
/* Name: area::checkout()
 * Function: Display the width, height, the first and the last element
 *	of p*/
/************************************************************************/
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

/************************************************************************/
/* Name: area::save2file(string name).
 * Function: Save the data of p in a file.
 * Description: name is the name of the file to be written.
 *	Data is written in rows. The first row is the highest row of p,
 *	i.e. p[height][]. Data of p of every time are separated by blank line.*/
/************************************************************************/
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
