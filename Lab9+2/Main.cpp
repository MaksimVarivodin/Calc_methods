#define _USE_MATH_DEFINES 
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
using namespace std;
double fu(const double & x) {
	return sin(exp(x / 3) + x);
}

double h_rec_tr(const double & a, const double & b, const double & N) {
	return (b - a) / N;
}
double left_rec(const double & N, const double & a = 3.0, const  double &b = 8.0) {
	double h = h_rec_tr(a, b, N);
	double integral = 0;
	for (double i = 0; i <= N - 1; i++) {
		integral += fu(a + i * h);
	}
	return integral * h;
}
double mid_rec(const double & N, const double & a = 3.0, const  double &b = 8.0) {
	double h = h_rec_tr(a, b, N);
	double integral = 0;
	for (double i = 0; i <= N - 1; i++) {
		integral += fu(a + h*(i + 0.5));
	}
	return integral * h;
}
double right_rec(const double & N, const double & a = 3.0, const  double &b = 8.0) {
	double h = h_rec_tr(a, b, N);
	double integral = 0;
	for (double i = 1; i <= N; i++) {
		integral += fu(a + i * h);
	}
	return integral * h;
}
double R_rec(const double & a, const double & b, const double & N) {
	double h = h_rec_tr(a, b, N), M = 5.6936;
	h *= h;
	return (h* M * (b - a)) / 24;
}
double trap(const double & N, const double & a = 3.0, const  double &b = 8.0) {
	double h = h_rec_tr(a, b, N);
	int res = N / 2.0;
	if(N / 2.0 - res > 0) return 0.0;
	double integral = (fu(a)+ fu(b)) / 2.0;
	for (double i = 1; i <= N - 1; i++) {
		integral += fu(a + i * h)*(a + (i + 1) * h - (a + (i - 1) * h) );
	}
	return integral* h;
}
double R_trap(const double & a, const double & b, const double & N) {
	double h = h_rec_tr(a, b, N), M = 28;
	h *= h;
	return (h*M*(b - a)) / 12;
}
double h_parabola(const double & a, const double & b, const double & N) {
	return (b - a) / (2.0 * N);
}
double parab(const double & N, const double & a = 3.0, const  double &b = 8.0) {
	double h = h_parabola(a, b, N);
	double buf1 = 0, buf2 = 0;
	for (double i = 1; i <= N; i++) {
		buf1 += 4.0 * fu(a + h * (2.0 * i - 1.0));
		if(i< N) buf2 += 2.0 * fu(a + 2.0 * i * h);
	}
	return (buf1+ buf2 + fu(a) + fu(b)) * (h / 3.0);
}
double R_parab(const double & a, const double & b, const double & N) {
	return (h_parabola(a, b, N)* 850.5921*(b - a)) / 180;
}
vector <double>
g_4_x
{
-0.8611363,
-0.3399810,
0.3399810,
0.8611363
},
g_4_c
{
0.3478548,
0.6521452,
0.6521452,
0.3478548
},
g_5_x
{ 
-0.9061798,
- 0.5384693,
0,
0.5384693,
0.5384693
},
g_5_c
{ 
0.2369269,
0.4786287,
0.5688888,
0.4786287,
0.2369269
},
g_6_x
{
-0.9324695,
- 0.6612093,
- 0.2386192,
0.2386192,
0.6612093,
0.9324695 },
g_6_c
{
0.1713245,
0.3607616,
0.4679131,
0.4679131,
0.3607616,
0.1713245
},
c_4_x
{ 
-0.794655,
- 0.187593,
0.187593,
0.794655
},
c_5_x
{
-0.832498,
- 0.374541,
0,
0.374541,
0.832498
},
c_6_x
{
-0.866247,
- 0.422519,
- 0.266636,
0.266636,
0.422519,
0.866247
},
k_4_h
{
	7.0 / 45.0,
	32.0 / 45.0,
	12.0 / 45.0, 
	32.0 / 45.0,
	7.0 / 45.0
},
k_5_h
{
	19.0 / 144.0,
	75.0 / 144.0,
	50.0 / 144.0,
	50.0 / 144.0,
	75.0 / 144.0,
	19.0 / 144.0
};

double gauss(const double & N, const double & a = 3.0, const  double &b = 8.0) {
	double res = 0, x = 0, c = 0;
	for (double i = 1; i <= N; i++) {
		if (N == 4) {
			x = g_4_x[i - 1];
			c = g_4_c[i - 1];
		}
		if (N == 5) {
			x = g_5_x[i - 1];
			c = g_5_c[i - 1];
		}
		if (N == 6) {
			x = g_6_x[i - 1];
			c = g_6_c[i - 1];
		}
		double buf1 = b * (1.0 + x), buf2 = a * (1.0 - x);
		res += c * fu((buf1 + buf2)/2.0);
	}
	return (res * (b - a)) / 2.0;
}
double chebyshev(const double & N, const double & a = 3.0, const  double &b = 8.0) {
	double res = 0;
	for (double i = 1; i <= N; i++) {
		if (N == 4)
			res += fu((c_4_x[i - 1] * (b - a) + a + b) /2.0);
		if (N == 5)
			res += fu((c_5_x[i - 1] * (b - a) + a + b) / 2.0);
		if (N == 6)
			res += fu((c_6_x[i - 1] * (b - a) + a + b) / 2.0);
		
	}
	
	return (res * (b - a) )/ N;
}
double kotesa(const double & N, const double & a = 3.0, const  double &b = 8.0) {
	double res = 0, h = (b - a)/N;

	for (double i = 0; i <= N; i++) {
		if (N == 4)
			res += k_4_h[i] * fu(a+ (i* h));
		if (N == 5)
			res += k_5_h[i] * fu(a + (i* h));
		
	}
	return (res * (b - a)) / 2.0;
}
int main()
{
	vector <string> names{
	"l_r mtd:",
	"m_r mtd:",
	"r_r mtd:",
	"trp mtd:",
	"prb mtd:",
	"gss mtd:",
	"chb mtd:", 
	"kts mtd:"
	};

	for (string var : names) {
		cout << setw(var.length() + 4) << var;
	}
	cout << endl;
	for (double i = 4; i < 7; i++)
	{
		vector < double > vec{
		left_rec(i),
		mid_rec(i),
		right_rec(i),
		trap(i),
		parab(i),
		gauss(i),
		chebyshev(i)
		};
		cout << "i: " << i;
		for (double var : vec) {
			cout << setw(to_string(var).length() + 3) << var;
		}
		if (i < 6){
			double buf = kotesa(i);
			cout << setw(to_string(buf).length() + 3) << buf;
		}
		cout << endl;
	}
	
		
}