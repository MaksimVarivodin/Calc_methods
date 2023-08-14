#pragma once
#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

class Matrix {
	double determinant;
	Matrix L_calc();
	Matrix U_calc();
	double det_calc();
	double(Matrix::*double_fu)(const int &, const int &, int &);// указатель на функции типа флоат
	void norm_assert(const int & Xiter, const int & Yiter, const bool & colorlin);
	// colorlin = column or line
	// Xiter = start iterator for a column
	// Yiter = start iterator for a line
	// в зависимости от того какой выбран режим(строка или колонка), \
	   мы можем по разному просматривать матрицу и ее вектора
	double mode_switcher(const int & Xiter, const int &Yîter, const bool & colorlin);// переключатель режимов (колонка или строка)
	double N_v_m_inside(const int & Xiter, const int &Yiter, int & size);// поиск максимального
	double N_v_l_inside(const int & Xiter, const int &Yiter, int & size);// сумма по модулю
	double N_v_e_inside(const int & Xiter, const int &Yiter, int & size);// Эвклидова норма 
	int sizeX;
	int sizeY;
	vector<double> mtrx;
	void line_swap(const int &i = 0, const int & j = 0);
	bool converge(Matrix xk, Matrix xkp, double eps);
	double okr(double x, double eps);
	bool diagonal(Matrix a);
	// новые
public:
	Matrix()
	{
		sizeX = 0;
		sizeY = 0;
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
	Matrix(vector<double>v, const int & sizeX, const int & sizeY) {
		mtrx = v;
		this->sizeX = sizeX;
		this->sizeY = sizeY;
	}
	int getX() { return sizeX; }
	int getY() { return sizeY; }
	vector<double>& get_mtrx() { return mtrx; }
	void Enter_check(int & a, const int & max = 1000, const int & min = -1000);
	Matrix& operator=(const Matrix& other);
	friend istream & operator>> (istream & in, Matrix & m);
	friend ostream & operator<< (ostream & out, const Matrix & m);
	Matrix operator^(const int & a);// возведение
	Matrix operator*(const Matrix& other);// умножение матриц
	Matrix operator*(const double& skalyar);// умножение матрицы на скаляр
	Matrix sum(const Matrix& other);// сумма двух матриц
	Matrix T();// транспонированная матрица
	Matrix E(const int & poryadok);// еденичная матрица заданного порядка
	double tr();// след
	double N_v_m_outside(const bool & colorlin = false, const int & Xiter = 0, const int &Yiter = 0);// поиск максимального
	double N_v_l_outside(const bool & colorlin = false, const int & Xiter = 0, const int &Yiter = 0);// сумма по модулю
	double N_v_e_outside(const bool & colorlin = false, const int & Xiter = 0, const int &Yiter = 0);// Эвклидова норма 
	double N_m_m();// манхеттенская норма
	double N_m_c();// норма чебышева
	double N_m_e();// евклидова норма
	double N_m_s();// спектральная норма
	// новые методы

	Matrix get_L() {
		return L_calc();
	}
	Matrix get_U() {
		return U_calc();
	}
	double get_det() {
		return determinant = det_calc();
	}
	// новые
	double cond1();
	double cond2();
	double condI();
	double condE();
	Matrix m_minus_1();// обратная матрица
	Matrix X(Matrix d);


	// решение слау при помощи метода Гаусса
	Matrix Zejdel(Matrix b);
};
