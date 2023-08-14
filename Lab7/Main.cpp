#include "Matrix_cl.h"
int main() {
	vector<float > a{2, -2,
		          -1, 3, 
		           1, 1,
		           2, -4}, 
	           	b{2, -4, 3, -4, 5, 
		         -1, -2, 3, 5, -4}, 
	         	c{2, 1, -4, 1, 4,
	              1, 2, 0, 2, 1, 
                  -1, 0, -6, 0, -1, 
	              1, 2, 0, 4, 1, 
	              4, 1, -2, 1, 4}, 
		        d{0, 1, -4, 1, 2, 
	              1, 4, -2, 2, 1,
	              2, 0, 6, 2, -1, 
	              1, 2, -2, 4, 1, 
	              4, 1, -2, 1, 4};
	int a_sizeY = 4, a_sizeX = 2,
		b_sizeY = 2, b_sizeX = 5,
		c_sizeY = 5, c_sizeX = 5,
		d_sizeY = 5, d_sizeX = 5;
	Matrix A{a, a_sizeX, a_sizeY};
	Matrix B{b, b_sizeX, b_sizeY};
	Matrix C{c, c_sizeX, c_sizeY};
	Matrix D{d, d_sizeX, d_sizeY};
	cout << "A^(T) " << A.T() << endl;
	cout << "B^(T) " << B.T()<< endl;
	cout << "C^(T) " << C.T()<< endl;
	cout << "D^(T) " << D.T()<< endl;
	cout << "A x B " << (A* B) << endl;
	cout << "B x A " << (B* A) << endl;
	cout << "C x D " << (C* D) << endl;
	cout << "D x C " << (D* C) << endl;
	cout << "A + B " << (A.sum( B)) << endl;
	cout << "C + D " << (C.sum( D)) << endl;
	cout << "C^(5) " << (C^(5)) << endl;
	cout << "D^(3) " << (D^(3)) << endl;
	cout << "A * 0.125 " << (A*0.125) << endl;
	cout << "B * -2 " << (B*(-2)) << endl;
	cout << "tr(B) " << (B.tr()) << endl;
	cout << "tr(D) " << (D.tr()) << endl;

}