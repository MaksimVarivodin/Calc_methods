#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;
class NEWT {
	clock_t start;
	clock_t end;
	double* x;
	int size;
	double *y;
	double koord_x;
	double koord_y;
public:
	NEWT()
	{
		this->x = NULL;
		this->y = NULL;
		this->koord_x = 0;
		this->koord_y = 0;
		this->size = 0;
		
		this->start = 0;
		this->end = 0;
	}
	NEWT(const int & size0):NEWT()
	{
		this->size = size0;	
		this->x = new double [size];
		this->y = new double [size];
	}
	~NEWT() {
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
	friend ostream& operator<< (ostream & out, NEWT& a );
	friend ofstream& operator<< (ofstream & out, NEWT& a);
	friend ifstream& operator>> (ifstream & in, NEWT& a);
	friend istream& operator>> (istream & in, NEWT& a);
	double fraction(const int & pos, const double & x0 = 0);
	void polynom();
	void put_m_info_to_file();
};