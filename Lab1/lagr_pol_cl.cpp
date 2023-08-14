#include "lagr_pol_cl.h"

double NEWT::fraction(const int & pos, const double & x0)
{
	double M = 1;	
	int buf = x0;
	for (int i = 0; i < size; i++)
	{
		if (i == pos) { i++; }
		switch (buf) {
		case 0:// ��� �����������
			M *= this->x[pos] - x[i];
			break;
		default:// ��� ���������
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

	cerr << "������� ������� ��������:" << endl;
	cin >> b;
	if (cin.fail()) return out;
	b++;
	a.x = new double[b];
	a.y = new double[b];
	a.size = b;
	for (int i = 0; i < a.size; i++) {
		cerr << "������� x" << i << " y"<< i<< " : "<< endl;		
		cin >> a.x[i];
		if (cin.fail()) return out;
		cin >> a.y[i];
		if (cin.fail()) return out;
	}
	cerr << "������� ���������� � ��� ������ �: " << endl;
	cin >> a.koord_x;

	return out;
}
ostream& operator<< (ostream & out, NEWT& a) {
	cerr << "������� ��������: "<< a.size - 1<< endl;
	cerr << "���������� �: " << a.koord_x << endl;
	cerr << "���������� �: " << a.koord_y << endl;
	for (int i = 0; i < a.size; i++)
	{
		cout << "��������� x" << i << " y" << i << " : "<< endl;
		cout << a.x[i]<< endl;
		cout <<  a.y[i]<< endl;
	}
	double buf = double(a.end - a.start);

	cout << "����� �������: ���: " << buf << endl;
	return out;

}