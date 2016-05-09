#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

/************************************************************************/
/*Name: area
 *Description: This class is the base for a 2D array.
 *member: 
 *	width: the width of the 2D array.
 *	height: the height of the 2D array.
 *	grid_num: the number of elements of p.
 *	p: store the data of the 2D array.
 *	filename: the name of the file store the data of the 2D array
 *	PI, mu, epsilon: some constants.
 *	area(int, int, string): initialize every member in this class.
 *	checkout(): display essential info about this class.
 *	save2file(string): save info to whose name is equal to filename.*/
/************************************************************************/
class area
{
public:
	int width;
	int height;
	int grid_num;
	vector<float> p;
	vector<float>::iterator itb, ite;
	string filename;
	const float PI = 3.14159265f;
	const float mu = (4.0 * PI) * 1e-7f;
	const float epsilon = 8.85e-12f;

public:
	area(int, int, string);
	void checkout();
	void save2file(string);
	~area();
};