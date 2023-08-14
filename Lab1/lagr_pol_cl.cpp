#include "lagr_pol_cl.h"

double NEWT::fraction(const int & pos, const double & x0)
{
	double M = 1;	
	int buf = x0;
	for (int i = 0; i < size; i++)
	{
		if (i == pos) { i++; }
		switch (buf) {
		case 0:// для знаменателя
			M *= this->x[pos] - x[i];
			break;
		default:// для числителя
			M *= x0 - x[i];
			break;
		}
	}
	return  M;
}

void NEWT::polynom()
{
	double res = 0;
	this->start = clock();
	for (int i = 0;i< this->size;i++) {
		double buf = this->y[i];
		res += ( buf* fraction(i, this->koord_x))/ fraction(i);
	}
	this->end = clock();
	this->koord_y = res;
}

void NEWT::put_m_info_to_file()
{

}

istream & operator>>(istream & out, NEWT & a)
{
	// TODO: insert return statement here
	int b;

	cerr << "Введите степень полинома:" << endl;
	cin >> b;
	if (cin.fail()) return out;
	b++;
	a.x = new double[b];
	a.y = new double[b];
	a.size = b;
	for (int i = 0; i < a.size; i++) {
		cerr << "Введите x" << i << " y"<< i<< " : "<< endl;		
		cin >> a.x[i];
		if (cin.fail()) return out;
		cin >> a.y[i];
		if (cin.fail()) return out;
	}
	cerr << "Введите координату х для поиска у: " << endl;
	cin >> a.koord_x;

	return out;
}
ostream& operator<< (ostream & out, NEWT& a) {
	cerr << "Степень полинома: "<< a.size - 1<< endl;
	cerr << "Координата х: " << a.koord_x << endl;
	cerr << "Координата у: " << a.koord_y << endl;
	for (int i = 0; i < a.size; i++)
	{
		cout << "Введенные x" << i << " y" << i << " : "<< endl;
		cout << a.x[i]<< endl;
		cout <<  a.y[i]<< endl;
	}
	double buf = double(a.end - a.start);

	cout << "Время расчета: сек: " << buf << endl;
	return out;

}