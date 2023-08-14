#include "Matrix_cl.h"

void Matrix::Enter_check(int & a, const int & max, const int & min)
{
	bool ready = false;
	string er = "Число не входит в диапазон значений!\n Введите число еще раз: ";
	while (!ready) // цикл продолжается до тех пор, пока пользователь не введет корректное значение
	{

		string b;
		// функция для ввода строки с пробелами

		getline(cin, b);

		try
		{
			a = stoi(b);
			if (a < min || a> max)  throw er;
			ready = true;
		}
		catch (const string & er)
		{
			cout << er;
		}
		catch (const std::exception&)
		{
			cout << "Замечены сторонние символы!\n Введите число еще раз: ";
		}
	}
}

Matrix & Matrix::operator=(const Matrix & other)
{
	mtrx = other.mtrx;
	sizeX = other.sizeX;
	sizeY = other.sizeY;
	return *this;
}

Matrix Matrix::operator^(const int & a)
{
	if (a > 0) {	
		Matrix buf{ *this };
		for (int i = 0; i < a; i++) {
			buf = buf * a;
		}
		return buf;
	}
	else if (a == 0) { 
		Matrix b{sizeX, sizeY};
		b = b.E(b.sizeX);		
		return b;
	}
	else {
		Matrix buf{ 0, 0 };
		return buf;
	}
	
}

Matrix  Matrix::operator*(const Matrix & other)
{

	Matrix a{ other.sizeX, other.sizeY };
	if (sizeX == other.sizeY && sizeY == other.sizeX) {

		// TODO: вставьте здесь оператор return
		for (int i = 0; i < a.sizeY; i++) {
			for (int j = 0; j < a.sizeX; j++) {
				for (int k = 0; k < sizeX; k++)
				{
					a.mtrx[i*a.sizeX + j] += mtrx[i*a.sizeX + k] * other.mtrx[k*other.sizeX + j];
				}
			}
		}
		return a;
	}
	else if (sizeX == other.sizeX && other.sizeY == 1) {
		Matrix buf{ other.mtrx, other.sizeY, other.sizeX };
		for (int i = 0; i < sizeY; i++) {
			for (int j = 0; j < sizeX; j++) {
				a.mtrx[i] += mtrx[i* sizeX + j] * buf.mtrx[j];
			}
		}
		return a;
	}
	else {
		cout << "Different sizes";
		Matrix buf{ 0, 0 };
		return buf;
	}

}

Matrix  Matrix::operator*(const float & skalyar)
{
	for (int i = 0; i < sizeY; i++)
		for (int j = 0; j < sizeX; j++)
			mtrx[i* sizeX + j] *= skalyar;
	return *this;
}

Matrix Matrix::sum(const Matrix & other)
{
	if (sizeX == other.sizeX && sizeY == other.sizeY) {
		Matrix A{ sizeX, sizeY };
		for (int i = 0; i < sizeY; i++) {
			int j = 0;
			for (; j < sizeX; j++) {
					A.mtrx[i*sizeX + j] = 0;
					A.mtrx[i*sizeX + j] = mtrx[i*sizeX + j] + other.mtrx[i*other.sizeX + j];
				}

			
		}
		return A;
	}
	else {
		Matrix buf{ 0, 0 };
		return buf;
	}
}



Matrix Matrix::T()
{
	Matrix a{ this->sizeY, this->sizeX };
	for (int i = 0; i < sizeY; i++) {
		for (int j = 0; j < sizeX; j++)
			a.mtrx[j*sizeY + i] = mtrx[i*sizeX + j];
	}
	return a;
}

Matrix Matrix::E(const int & poryadok)
{
	Matrix a{poryadok, poryadok};
	a.mtrx.resize(poryadok* poryadok);
	for (int i = 0; i < poryadok; i++) {
		for (int j = 0; j < poryadok; j++) {
			if (i == j)a.mtrx[i*poryadok + j] = 1;
			else 
			a.mtrx[i*poryadok+j] = 0;

		}
	}		
	return a;
}

float Matrix::tr()
{
	int i = 0, j = 0;
	
	if (sizeY == sizeX) {
	    float  to_return = 0;
		for (; i < sizeY; i++) {
			j= i;
			to_return +=mtrx[i*sizeX + j];

		}
		return to_return;
	}
	else return 0;
}

istream & operator>>(istream & in, Matrix & m)
{
	cout << "Enter data of the matrix: " << endl;
	for(int i  = 0; i< m.sizeY; i++)
		for (int j = 0; j < m.sizeX; j++)
		{
			cout << "Enter i:"<< i<< " j: " << j;
			int buf = 0;
			m.Enter_check(buf);
			m.mtrx[i * m.sizeX + j] = buf;
		}
	return in;
}

ostream & operator<<(ostream & out,const  Matrix & m)
{
	cout << "Matrix: " << endl;
	for (int i = 0; i < m.sizeY; i++) {
		for (int j = 0; j < m.sizeX; j++)
		{
			cout << m.mtrx[i * m.sizeX + j]<< "\t";			
		}
		cout << endl;
	}
		
	return out;
}

