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
public:
	area(int, int, string);
	void checkout();
	void save2file(string);
	//virtual void real_pos(int, int);
	~area();
};