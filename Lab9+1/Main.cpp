#include "Matrix4_cl.h"
int main() {

	vector <double>
	
		a{
		-3, 9, 3, 6, 
		-5, 8, 2, 7, 
		4, -5, -3, -2, 
		7, -8, -4, -5
		},
		b{
		 3, 2, 2, 2, 
	     9, -8, 5, 10, 
	     5, -8, 5, 8, 
	     6, -5, 4, 7},
		c{ 1.0, 1.0 / 4.0, -1.0 / 4.0, 0.1,
		  -0.1, 0.9, 0.2, 1.0  / 6.0,
		  -0.2, -0.1/ 6.0, 0.8, 0.2,
		  0.2, 0.1, -0.1, 1.0},
		d{ 1.1, 7.0/6.0, 38.0 / 60.0, 1.2 },
		x{ 1.1, 1.4, 1.9, 1.1 };
	int a_sizeY = 4, a_sizeX = 4,
		b_sizeY = 4, b_sizeX = 4,
		c_sizeY = 4, c_sizeX = 4,
		d_sizeY = 1, d_sizeX = 4,
		x_sizeY = 1, x_sizeX = 4;
	Matrix A{ a, a_sizeX, a_sizeY };
	Matrix B{ b, b_sizeX, b_sizeY };
	Matrix C{ c, c_sizeX, c_sizeY };
	Matrix D{ d, d_sizeX, d_sizeY };
	Matrix X{ x, x_sizeX, x_sizeY };
	double cond1_A = A.cond1(),
		cond1_B = B.cond1(),
		m1_A = 0, m1_B = 0,
		m2_A = 0, m2_B = 0,
		m3_A = 0, m3_B = 0;

	cout << "X for A "<< (A.X(X))<< endl;
	cout << "cond1_A " << cond1_A << endl;
	cout<< "X for B "<< (B.X(X)) << endl;
	cout << "cond1_B " << cond1_B << endl;

	A.get_mtrx()[0] += 0.1;
	B.get_mtrx()[0] += 0.1;
	m1_A = (cond1_A* 0.1) / A.N_m_m();
	m1_B = (cond1_B* 0.1) / B.N_m_m();
	cout << "X for A + 0.1 " << A.X(X) << endl;
	cout << "modificaction 1 condition "<< m1_A << endl;
	cout << "X for B + 0.1 " << B.X(X) << endl;
	cout << "modificaction 1 condition " << m1_B << endl;


	A.get_mtrx() =a;
	B.get_mtrx() =b;
	X.get_mtrx()[0] += 0.1;
	m2_A = (cond1_A* 0.1) / X.N_v_m_outside();
	m2_B = (cond1_B* 0.1) / X.N_v_m_outside();
	cout << "X for A,  X + 0.1 " << A.X(X) << endl;
	cout << "modificaction 2 condition " << m1_A << endl;
	cout << "X for B,  X + 0.1 " << B.X(X) << endl;
	cout << "modificaction 2 condition " << m1_B << endl;


	A.get_mtrx()[0] += 0.1;
	B.get_mtrx()[0] += 0.1;
	X.get_mtrx()[0] += 0.1;
	m3_A = (cond1_A* 0.1* (A.N_m_m()+ X.N_v_m_outside())) / (A.N_m_m() * X.N_v_m_outside());
	m3_B = (cond1_B* 0.1* (B.N_m_m() + X.N_v_m_outside())) / (B.N_m_m() * X.N_v_m_outside());
	cout << "X for A + 0.1 , X + 0.1" << A.X(X) << endl;
	cout << "modificaction 3 condition " << m2_A << endl;
	cout << "X for B + 0.1 , X + 0.1" << B.X(X) << endl;
	cout << "modificaction 3 condition " << m2_B << endl;

	Matrix buf1 = C.X(D);
	Matrix buf2 = C.Zejdel(D);
	cout << buf1 << endl;
	cout << buf2 << endl;
	cout << "vector nevjazki for Prjamoj " << D.sum((C * buf1) * (-1.0))<< endl;
	cout << "vector nevjazki for Zejdel " << D.sum((C * buf2) * (-1.0)) << endl;
}