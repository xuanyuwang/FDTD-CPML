#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

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