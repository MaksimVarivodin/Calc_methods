#include "geom_pol.h"
int main() {
	const double pi = 3.14;
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	int size = 5;
	double *masX = new double [size];
	for (int i = 0; i < size; i++)
		masX[i] = i * 2 * pi / 5;
	double *masY = new double[size] {8.6, -12.2, -2.8, 2.7, -5.3};
	Geom arr[5]{ { masX ,masY, size },{ masX ,masY, size },{ masX ,masY, size },{ masX ,masY, size },{ masX ,masY, size } };
	for (int i = 0; i < size; i++) {
		double add = (1 / 5)*pi;
		arr[i].x_setter(masX[i] + add);
	}
	for (int i = 0; i < size; i++)
		arr[i].polynom();
	for (int i = 0; i < size; i++)
		cout << arr[i];
	system("pause");
	delete[] masX;
	delete[] masY;
}