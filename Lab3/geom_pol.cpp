#include "geom_pol.h"
void Geom::x_setter(const double & x0)
{
	koord_x = x0;
}
double Geom::fraction(const int & pos, const double & x0)
{
	double M = 1;
	int buf = x0;
	for (int i = 0; i < size; i++)
	{
		if (i == pos) { i++; }
		switch (buf) {
		case 0:// для знаменателя
			M *= sin((this->x[pos] - x[i]) / 2);
			break;
		default:// для числителя
			M *= sin(x0 - x[i] / 2);
			break;
		}
	}
	return  M;
}

void Geom::polynom()
{
	double res = 0;
	this->start = clock();
	for (int i = 0; i< this->size; i++) {
		double buf = this->y[i];
		res += (buf* fraction(i, this->koord_x)) / fraction(i);
	}
	this->end = clock();
	this->koord_y = res;
}


istream & operator>>(istream & out, Geom & a)
{

	cerr << "Введите координату х для поиска у: " << endl;
	cin >> a.koord_x;

	return out;
}
ostream& operator<< (ostream & out, Geom& a) {
	
	
	for (int i = 0; i < a.size; i++)
	{
		cout << "Введенные x" << i << " y" << i << " : " << endl;
		cout << a.x[i] << endl;
		cout << a.y[i] << endl;
	}
	cerr << "Координата х: " << a.koord_x << endl;
	cerr << "Координата у: " << a.koord_y << endl;
	
	return out;

}