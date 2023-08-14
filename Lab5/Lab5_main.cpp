#include <iostream>
#include <cmath>
#include <limits>
#include <Windows.h> 
using namespace std;
double fu(const double & x) {
	double var = x * x * x;
	var -= 2 * x * x;
	var -= 6 * x;
	var--;
	return var;
}
double March(const double & a) {
	return  3 * pow(a, 2) - 4 * a - 6;
}
double March2(const double& a) {
	return 6 * a - 4;
}
double Dehotomiya_method(const double & a, const double & b, double& accuracy) {
	double buf_a = a, buf_b = b, buf_c = 0, fu_a = 0, fu_b = 0, fu_c = 0;
	int iter = 0;
	while (buf_b - buf_a >= accuracy && iter < 1000) {		
		fu_a = fu(buf_a);
		fu_b = fu(buf_b);
		buf_c = (buf_a + buf_b) / 2;
		fu_c = fu(buf_c);
		if (fu_a * fu_c <= 0) {
			buf_b = buf_c;
		}
		else {
			buf_a = buf_c;
		}
		iter++;
	}
	buf_c = (buf_a + buf_b) / 2;
	cout << __FUNCTION__ << " count: " << iter << endl;
	cout << __FUNCTION__ << " A: " << a << endl;
	cout << __FUNCTION__ << " B: " << b << endl;
	cout << __FUNCTION__ << " C: " << buf_c << endl;
	return buf_c;
}

double Khord_method( const double &a ,const double &b,const double & accuracy)
{
	double fu1 = 0, fu2 = 0, 
		x0 = (a + b) / 2, x1 = 0;
	int iter = 0;
	
	fu1 = March(x0);
	fu2 = March2(x0);
	if (fu1 * fu2 < 0) {
		x1 = b;
		while(fabs(x1 - x0) >= accuracy){
			iter++;
			x0 = x1;
			x1 = x0 - (fu(x0)*(x0 - a)) / (fu(x0) - fu(a));
			
			cout << "fabs(x1 - x0): " << fabs(x1 - x0) << endl;
		}
	}
	else if (fu1 * fu2 > 0) {
		x1 = a;
		while (fabs(x1 - x0) >= accuracy) {
			iter++;
			x0 = x1;
			x1 = x0 - (fu(x0)*(b - x0)) / (fu(b) - fu(x0));
			
			cout << "fabs(x1 - x0): " << fabs(x1 - x0) << endl;
		}
	}
	cout << __FUNCTION__ << " count: " << iter<< endl;
	cout << __FUNCTION__ << " A: " << a << endl;
	cout << __FUNCTION__ << " B: " << b << endl;
	cout << __FUNCTION__ << " C: " << x1 << endl;
	return x1;
}
double Newton_method( const double & a, const double &b, const double &accuracy) {

	double x0 = 0, x1 = 0;
	int iter = 0;

	if (fu(a) * March2(a) > 0) {
		x1 = a;
	}
	else if (fu(b) * March2(b) > 0) {
		x1 = b;		
	}
	while (fabs(x1 - x0) >= accuracy || March(x1) == numeric_limits<double>::infinity()) {
		iter++;
		x0 = x1;
		x1 = x0 - (fu(x0)/March(x0));
		cout << "fabs(x1 - x0): " << fabs(x1 - x0) << endl;
	}
	cout << __FUNCTION__ << " count: " << iter << endl;
	cout << __FUNCTION__ << " A: " << a << endl;
	cout << __FUNCTION__ << " B: " << b << endl;
	cout << __FUNCTION__ << " C: " << x1 << endl;
	return x1;
	
}
double Secant_method(const double &a, const  double &b, const double & accuracy) {
	double x0 = 0, x1 = a, x2 = b;

	int iter = 0;
	while (fabs(x1 - x0) >= accuracy) {
		iter++;
		x0 = x1;
		x1 = x2;
		x2 = x1 - (fu(x1)*(x1 - x0)) / (fu(x1) - fu(x0));
		cout << "fabs(x1 - x0): " << fabs(x1 - x0) << endl;
	}
    cout << __FUNCTION__ << " count: " << iter << endl;
	cout << __FUNCTION__ << " A: " << a << endl;
	cout << __FUNCTION__ << " B: " << b << endl;
	cout << __FUNCTION__ << " C: " << x2 << endl;
	return x2;
}

double Combination_method(const double &a, const  double &b, const double & accuracy, const int & counter = 0) {
	
	double fu1 = 0, fu2 = 0,
		buf_a = a, buf_b = b;
	int iter = 0;

	fu1 = March((a+b)/2);
	fu2 = March2((a + b) / 2);
	if (fu1 * fu2 > 0) {
		while (fabs(buf_b - buf_a) >= accuracy) {
			iter++;
			buf_a = buf_a - (fu(buf_a)*(buf_b - buf_a)) / (fu(buf_b) - fu(buf_a));
			buf_b = buf_b - (fu(buf_b) / March(buf_b));
		}
	}
	else if (fu1 * fu2 < 0) {
		while (fabs(buf_b - buf_a) >= accuracy) {
			iter++;
			double buf_var = buf_a;
			buf_a = buf_a - (fu(buf_a) / March(buf_a));
			buf_b = buf_b - (fu(buf_b)*(buf_b - buf_var)) / (fu(buf_b) - fu(buf_var));			
		}
	}
	cout << __FUNCTION__ << " count: " << iter << endl;
	cout << __FUNCTION__ << " A: " << a << endl;
	cout << __FUNCTION__ << " B: " << b << endl;
	cout << __FUNCTION__ << " C: " << (buf_a + buf_b)/2 << endl;
	return (buf_a + buf_b) / 2;
}
double gi_ot_x(const double & x) {
	return  1 / (x* x - 2 * x - 6);
}
double Iter_method(const double &a, const  double &b, const double & accuracy) {

	double x0 = (a + b)/2, x1 = gi_ot_x(x0);
	int iter = 0;
	while (abs(x1 - x0) > accuracy) {
		iter++;
		cout << "x0 " << x0 << endl;
		x0 = x1;
		if (gi_ot_x(x0) == numeric_limits<double>::infinity())
			break;
		x1 = gi_ot_x(x0);
	}
	cout << __FUNCTION__ << " count: " << iter << endl;
	cout << __FUNCTION__ << " A: " << a << endl;
	cout << __FUNCTION__ << " B: " << b << endl;
	cout << __FUNCTION__ << " C: " << x1 << endl;
	return x1;
}
int main() {
	cout << "There are three roots in intervals:" << endl;
	cout << "-2.5; -1" << endl;
	cout << "-0.8; 2" << endl;
	cout << "3; 4" << endl;
	
	double a = 3, b = 4,
		   c = -0.5, d = 2, 
		   e  = -2.5 , f  = - 1,
		   ac = 0.001;

    Dehotomiya_method( c, d, ac);
	Khord_method(c, d, ac);
	Newton_method(c, d, ac);
	Secant_method(c,d, ac);
	Combination_method(c, d, ac);
	Iter_method(c,d, ac);
	
	system("pause");
}
