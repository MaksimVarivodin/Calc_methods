#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;
class Geom {
	clock_t start;
	clock_t end;
	double* x;
	int size;
	double *y;
	double koord_x;
	double koord_y;
public:
	Geom()
	{
		this->x = NULL;
		this->y = NULL;
		this->koord_x = 0;
		this->koord_y = 0;
		this->size = 0;
		
		this->start = 0;
		this->end = 0;
	}
	Geom(const double* arrX,const double* arrY, const int & size0)
	{
		this->size = size0;
		x = new double[size];
		y = new double[size];
		for (int i = 0; i < size; i++)
		{
			x[i] = arrX[i];
			y[i] = arrY[i];
		}
	}
	Geom(const int & size0):Geom()
	{
		this->size = size0;	
		this->x = new double [size];
		this->y = new double [size];
	}
	~Geom() {
		if (size > 0)
		{
			delete[] x;
			delete[] y;
		}
		
		this->koord_x = 0;
		this->koord_y = 0;
		this->size = 0;		
		this->start = 0;
		this->end = 0;
	}
	void x_setter(const double &x0);
	friend ostream& operator<< (ostream & out, Geom& a );
	friend ofstream& operator<< (ofstream & out, Geom& a);
	friend ifstream& operator>> (ifstream & in, Geom& a);
	friend istream& operator>> (istream & in, Geom& a);
	double fraction(const int & pos, const double & x0 = 0);
	void polynom();
};
