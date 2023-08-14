#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;
 

class Matrix {
protected:
	int sizeX;
	int sizeY;
	vector<float> mtrx;
public:
	Matrix()
	{
		cout << "Enter sizeX: "; Enter_check(sizeX, 100, 0);
		cout << "Enter sizeY: "; Enter_check(sizeY, 100, 0);
		mtrx.resize(sizeX * sizeY);
	}
	Matrix(const int & a, const int& b)
	{
		sizeX = a;
		sizeY = b;
		mtrx.resize(a*b);
	}
	Matrix(const Matrix &other)
	{
		*this = other;
	}
	Matrix(vector<float>v, const int & sizeX, const int & sizeY) {
		mtrx = v;
		this->sizeX = sizeX;
		this->sizeY = sizeY;
	}
	int getX(){ return sizeX; }
	int getY() { return sizeY; }
	vector<float>& get_mtrx() { return mtrx; }
	void Enter_check(int & a, const int & max  = 1000, const int & min = -1000);
	Matrix& operator=(const Matrix& other);
	friend istream & operator>> (istream & in, Matrix & m);
	friend ostream & operator<< (ostream & out, const Matrix & m);
	Matrix operator^(const int & a);// возведение
	Matrix operator*(const Matrix& other);// умножение матриц
	Matrix operator*(const float& skalyar);// умножение матрицы на скаляр
	Matrix sum(const Matrix& other);// сумма двух матриц
	Matrix T();// транспонированная матрица
	Matrix E(const int & poryadok);// еденичная матрица заданного порядка
	float tr();// след
};