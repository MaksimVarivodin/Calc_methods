#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
using namespace std;

double f1_x(const double & x, const double & y) {
	return 2.0 * x + y;
}
double f1_y(const double & x, const double & y) {
	return x - 4.0;
}
double f2_x(const double & x, const double & y) {
	return (1.0 / y) - 1.0;
}
double f2_y(const double & x, const double & y) {
	return (((-1)*x) / (y * y)) - 1.0;
}
double znam_stisk_oper(const double & x, const double & y) {
	return f2_x(x, y) * f1_y(x, y) - f2_y(x, y) * f1_x(x, y);
}
double alfa(const double & x, const double & y) {
	return f2_y(x, y) / znam_stisk_oper(x, y);
}
double beta(const double & x, const double & y) {
	return (-1.0 )*(f1_y(x, y) / znam_stisk_oper(x, y));
}
double gamma(const double & x, const double& y) {
	return (-1.0)*(f2_x(x, y) / znam_stisk_oper(x, y));
}
double delta(const double & x, const double & y) {
	return f1_x(x, y) / znam_stisk_oper(x, y);
}
void g_vec(vector<double> &vec) {
	double x0 = vec[0], 
		   y0 = vec[1], 
		   x = vec[2], 
		   y = vec[3], 
		   a = alfa(x0, y0), 
		   b = beta(x0, y0),
		   g = gamma(x0, y0),
		   d = delta(x0, y0);
	vec[0] = vec[2];
	vec[1] = vec[3];
	vec[2] = x + a * (x*x + y * x - 4 * y) + b * ((x / y) - x - y); 
    vec[3] = y + g * (x*x + y * x - 4 * y) + d * ((x / y) - x - y);
}
vector<double> iter(const double & x, const double & y, const double & E = 0.00001) {
	vector<double>vec1{ x, y, x, y };
	while (abs(vec1[0] - vec1[2]) < E ||  abs(vec1[1] - vec1[3]) < E) {
		g_vec(vec1);
	}
	return vector<double>{vec1[2], vec1[3]};
}
void ober_matr(double a[2][2])
{
	double det, aa;
	det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
	aa = a[0][0];
	a[0][0] = a[1][1] / det;
	a[1][1] = aa / det;
	aa = a[0][1];
	a[0][1] = -a[0][1] / det;
	a[1][0] = -a[1][0] / det;
}

vector < double > newton(double x, double y, const double & E = 0.00001)
{
	int i = 1;
	double a[2][2], dx, dy, b[2], norm;
	do
	{
		a[0][0] = f1_x(x, y);
		a[0][1] = f1_y(x, y);
		a[1][0] = f2_x(x, y);
		a[1][1] = f2_y(x, y);
		ober_matr(a);
		dx = -a[0][0] * (x*x + y * x - 4 * y) + -a[0][1] * ((x / y) - x - y);
		dy = -a[1][0] * (x*x + y * x - 4 * y) + -a[1][1] * ((x / y) - x - y);
		x = x + dx;
		y = y + dy;
		b[0] = (x*x + y * x - 4 * y);
		b[1] = ((x / y) - x - y);
		norm = sqrt(b[0] * b[0] + b[1] * b[1]);
		i++;
	} while (norm >= E);
	return vector <double> {x, y};
}
int main() {
	vector <string> names{ "Simple iteration method,", "Simple iteration method," , "Newton method,", "Newton method,", " root " , "x: ", "y: "};
	vector <double> iterat[4]{ iter(-4.5, 2.5), iter(1.2, 0.5), newton(-4.5, 2.5), newton(1.2, 0.5) };
	int i = 0, root = 1, it = 1;

	for (vector < double > var : iterat) {
		cout << setw(names[i].length()+ 10) << names[i] << names[4]<< root << endl;
		for (double var1 : var) {
			if (it % 2 == 0)cout << names[6];
			else cout << names[5];
			cout << var1 << endl;
			it++;
		}
		if (i == 1)root = 1;
		else root++;
		i++;
	}
	
	
}