#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

class area
{
public:
	int length;
	int grid_num;
	vector<float> p;
	vector<float>::iterator itb, ite;
	string filename;
public:
	area(int, string);
	void checkout();
	void save2file(string);
	~area();
};