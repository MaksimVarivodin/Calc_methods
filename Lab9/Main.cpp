#include "Matrix3_cl.h"
int main() {

	vector <double>
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
	int c_sizeY = 5, c_sizeX = 5,
		d_sizeY = 5, d_sizeX = 5;
	Matrix C{ c, c_sizeX, c_sizeY };
	Matrix D{ d, d_sizeX, d_sizeY };

	cout << "C "<< C << endl;
	cout << "L " << C.get_L() << endl;
	cout << "L x U "<< C.get_L() * C.get_U() << endl;
	cout << "Determinant " << C.get_det() << endl;
	cout << "D "<< D << endl;
	cout << "L " << D.get_L() << endl;
	cout << "Determinant " << D.get_det() << endl;

}