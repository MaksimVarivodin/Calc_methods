#include <iostream>
#include <cmath>
#include <Windows.h> 
using namespace std;
 // Параболическая апроксимация
double* S1
(double * Xi, double* Yi, const int & size) {
	double *s = new double[7]{
	0, 0, 0, 0, 0, 0, 0	
	};
	for (int i = 0; i < size; i++) {
		s[0] += Xi[i];
		s[1] += pow(Xi[i], 2);
		s[2] += pow(Xi[i], 3);
		s[3] += pow(Xi[i], 4);
		s[4] += Yi[i];
		s[5] += Xi[i] * Yi[i];
		s[6] += pow(Xi[i], 2) * Yi[i];
	}
	return s;
}
 
double determinant
   (const double & a1, const double &  a2, const double & a3,
	const double & a4, const double & a5, const double & a6, 
	const double & a7, const double & a8, const double & a9) {
	double sum1 = 0, sum2 = 0;
	sum1 = a1*a5* a9 + a4* a3* a8 + a7*a2*a6;
	sum2 = a3*a5* a7 + a4* a2* a9 + a8*a6*a1;
	return sum1 - sum2;
}
double delta0
(const  double &size, double *s) {
	return determinant(
		size, s[0], s[1],
		s[0], s[1], s[2],
		s[1], s[2], s[3]);
}

double delta1(double *s) {
	return determinant(
		s[4], s[0], s[1],
		s[5], s[1], s[2],
		s[6], s[2], s[3]);
}
double a0( double * s, const double & del) {
	return delta1(s) / del;
}
double delta2(const  double &size, double *s) {
	return determinant(
		size, s[4], s[1],
		s[0], s[5], s[2],
		s[1], s[6], s[3]);
}
double a1(const  double &size,double * s, const double & del) {
	return delta2(size,	s) / del;
}
double delta3(const  double &size, double *s) {
	return determinant(
		size, s[0], s[4],
		s[0], s[1], s[5],
		s[1], s[2], s[6]);
}
double a2
(	const  double &size, double * s, const double & del) {
	return delta3(size,	s ) /del;
}
double y_count
(	double * Xi, double *Yi, const double &size, const double &x)
{
	double *S = S1(Xi, Yi, size);	
	double del0 = delta0(size, S);
	double * arr = new double[3]{
	a0(S, del0), a1(size, S, del0),a2(size, S, del0)
	};
	double sum = arr[2] * pow(x, 2) + arr[1] * x + arr[0];
	delete[] S;
	delete[] arr;
	return sum;
}

// Линейная апроксимация
double a_count
(double * Xi, double* Yi, const int & size) 
{
		double a = 0;
	double sum1=0, sum2 = 0, sum3 = 0;
	for (int i = 0; i < size; i++) {
		sum1 += Xi[i] * Yi[i];
	}
	sum1*=size;
	
	for (int i = 0; i < size; i++) {
		sum2 += Xi[i];
		sum3 += Yi[i];
	}
	sum1 = sum1 - sum2*sum3;	
	sum3 =0;
	for (int i = 0; i < size; i++) {
		sum3 += pow(Xi[i], 2);
	}
	sum3 *= size;
	sum3 -= pow(sum2, 2);
	a = sum1 / sum3;
	return a;
}
double b_count
(const double &a, double * Xi, double * Yi, const int & size)
{
	double b = 0;
	double sum1 = 0, sum2 = 0, sum3 = 0;
	for (int i = 0; i < size; i++) {
		sum1 += (Yi[i] -a*Xi[i]);
	}
	b = sum1 / size;
	return b;
}



int main()
{
	int size = 10;
	double *Xi = new double[size]{ 7.5, 16.3, 0.2, 8.5, 16.8, 18.9, 7.7, 11.1, 9.6, 10.0 },
		*Yi = new double[size] {1.1, 6.4, 10.9, 12.5, 5.9, 15.0, 15.5, 2.1, 15.5, 16.8	};
	double x = 0;// не придумал как ввести
	double a = a_count(Xi, Yi, size);
	double b = b_count(a, Xi, Yi, size);
	double y = y_count(Xi, Yi, size, Xi[0]);
	cout << "A: " << a << endl << "B: " << b << endl;
	cout << y << endl;
	system("pause");
	delete [] Xi;
	delete[] Yi;
}