#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;
// y'= z
double f1(const double & z)
{
	return z;
}
// z' = (2 + 2/10)* z + (2/10)* y + e^(x)
double f2(const double & x, const double & y, double z)
{
	return 2.8 * z + 0.8 * y + exp(x);
}
void euler_original(
	 const double &xi, const double &yi, const double& zi,const double& xf,
	double& yf, double& zf)
{
	// y(k + 1) = y(k) + z * h 
	yf = yi + f1(zi) * (xf - xi);
	// z(k + 1) = z(k) + h *z'
	zf = zi + f2(xi, yi, zi) * (xf - xi);
}

void euler_2modific(
	double xi, double yi, double zi, double xf,
	double& yf, double& zf)
{ 
	// y* & z* 
	euler_original(xi, yi, zi, xf, yf, zf);
	// y(k + 1) = y(k) + (z(k)+ z(k+ 1))* h / 2
	yf = yi + (f1(zi) + f1(zf)) * 0.5 * (xf - xi);
	// z(k + 1) = z(k) + (z'(x, y, z)+ z'(x+ 1, y+ 1, z + 1))* h / 2
	zf = zi + (f2(xi, yi, zi) + f2(xf, yf, zf)) * 0.5 * (xf - xi);
}

void runge_kutt(
	double xi, double yi, double zi, double xf,
	double& yf, double& zf)
{
	double 
		h = xf - xi, 
		t = xi,
		k1x, k2x, k3x, k4x, k1v, k2v, k3v, k4v;


	k1x = h * f1(zi);
	k1v = h * f2(t, yi, zi);

	k2x = h * f1(zi + k1v / 2.0);
	k2v = h * f2(t + h / 2.0, yi + k1x / 2.0, zi + k1v / 2.0);

	k3x = h * f1( zi + k2v / 2.0);
	k3v = h * f2(t + h / 2.0, yi + k2x / 2.0, zi + k2v / 2.0);

	k4x = h * f1(zi + k3v);
	k4v = h * f2(t + h, yi + k3x, zi + k3v);

	yf = yi + (k1x + 2.0 * (k2x + k3x) + k4x) / 6.0;
	zf = zi + (k1v + 2.0 * (k2v + k3v) + k4v) / 6.0;
}
void print(int key, double xi, double yi, double zi, double hmax, double h) {

	double  xf, yf, zf;

	const string method[3] = { "Метод Эйлера","Метод Эйлера-Коши","Метод Рунге-Кутта 4-го порядка" };

	cout << setw(30) << method[key] << endl << endl;
	cout << setw(12) << "h" << setw(12) << "y" << setw(12) << "y'" << endl;
	cout << setw(12) << xi << setw(12) << yi << setw(12) << zi << endl;
	while (xi <= hmax)
	{
		xf = xi + h;
		switch (key) {
		case 0:euler_original(xi, yi, zi, xf, yf, zf);
			break;
		case 1:euler_2modific(xi, yi, zi, xf, yf, zf);
			break;
		case 2: runge_kutt(xi, yi, zi, xf, yf, zf);
			break;
		}

		cout << setw(12) << xf << setw(12) << yf << setw(12) << zf << endl;
		xi = xf;
		yi = yf;
		zi = zf;
	}
}



int main()
{
	setlocale(LC_ALL, "Russian");

	double xi, yi, zi, h, hmax;

	xi = 0.0;
	yi = 1.0;
	zi = 1.0;
	h = 0.1;
	hmax = 0.9;

	print(0, xi, yi, zi, hmax, h);
	cout << endl;
	print(1, xi, yi, zi, hmax, h);
	cout << endl;
	print(2, xi, yi, zi, hmax, h);

	system("pause");
	return 0;
}

