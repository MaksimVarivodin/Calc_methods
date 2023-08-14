#include "Matrix2_cl.h"
int main() {
	vector <double>
		x{1, 2, 3, 
	4, 5, 6, 
	7, 8, -8},
		a{ 2, -1, -5, 3, 4, 1, -2 },
		b{ -11, 2, 3, 5, 1 },
		c{ 2, 1, -4, 1, 4,
		  1, 2, 0, 2, 1,
		  -1, 0, -6, 0, -1,
		  1, 2, 0, 4, 1,
		  4, 1, -2, 1, 4 },
		d{ 0, 1, -4, 1, 2,
		  1, 4, -2, 2, 1,
		  2, 0, 6, 2, -1,
		  1, 2, -2, 4, 1,
		  4, 1, -2, 1, 4 };
	int x_sizeY = 3, x_sizeX = 3,
		a_sizeY = 1, a_sizeX = 7,
		b_sizeY = 5, b_sizeX = 1,
		c_sizeY = 5, c_sizeX = 5,
		d_sizeY = 5, d_sizeX = 5;
	Matrix X{ x, x_sizeX, x_sizeY };
	Matrix A{ a, a_sizeX, a_sizeY };
	Matrix B{ b, b_sizeX, b_sizeY };
	Matrix C{ c, c_sizeX, c_sizeY };
	Matrix D{ d, d_sizeX, d_sizeY };
	cout << "Max-norm: " << A.N_v_m_outside()<< " for vector A."<<endl;
	cout << "L-norm: " << A.N_v_l_outside() << " for vector A." << endl;
	cout << "Evklid norm: " << A.N_v_e_outside() << " for vector A." << endl;
	cout << endl;
	cout << "Max-norm: " << B.N_v_m_outside(true) << " for vector B." << endl;
	cout << "L-norm: " << B.N_v_l_outside(true) << " for vector B." << endl;
	cout << "Evklid norm: " << B.N_v_e_outside(true) << " for vector B." << endl;
	cout << endl;

	cout << "Manhattan norm: " << C.N_m_m() << " for matrix C." << endl;
	cout << "Chebyshev norm: " << C.N_m_c() << " for matrix C." << endl;
	cout << "Evklid-norm: " << C.N_m_e() << " for matrix C." << endl;
	cout << "Spectral-norm: " << C.N_m_s() << " for matrix C." << endl;
	cout << endl;
	cout << "Manhattan norm: " << D.N_m_m() << " for matrix D." << endl;
	cout << "Chebyshev norm: " << D.N_m_c() << " for matrix D." << endl;
	cout << "Evklid-norm: " << D.N_m_e() << " for matrix D." << endl;
	cout << "Spectral-norm: " << D.N_m_s() << " for matrix D." << endl;
	cout << endl;

	cout << "Manhattan norm: " << X.N_m_m() << " for matrix X." << endl;
	cout << "Chebyshev norm: " << X.N_m_c() << " for matrix X." << endl;
	cout << "Evklid-norm: " <<    X.N_m_e() << " for matrix X." << endl;
	cout << "Spectral-norm: " <<  X.N_m_s() << " for matrix X." << endl;
	cout << endl;
	// m4* m3* m2*m1 = M
} 