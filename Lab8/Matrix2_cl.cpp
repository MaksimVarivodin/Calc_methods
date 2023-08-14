#include "Matrix2_cl.h"

void Matrix::norm_assert (const int & Xiter, const int & Yiter, const bool & colorlin)
{
	bool correct_Y_and_X_for_line = (Xiter == 0 && (Yiter > -1 && Yiter < sizeY) && !colorlin), 
		 correct_Y_and_X_for_cols = (Yiter == 0 && (Xiter > -1 && Xiter < sizeX) && colorlin);
	assert(correct_Y_and_X_for_line || correct_Y_and_X_for_cols);
}
// если colorlin то столбцы, если нет то строки
double Matrix::mode_switcher(const int & Xiter, const int &Yiter, const bool & colorlin)
{
	int size = 1;
	if (colorlin)
		size = sizeX;
	return (this->*double_fu)(Xiter, Yiter, size);
}

double Matrix::N_v_m_outside(const bool & colorlin, const int & Xiter, const int & Yiter)
{
	norm_assert(Xiter, Yiter, colorlin);	
	this->double_fu = &Matrix::N_v_m_inside;
	return mode_switcher(Xiter, Yiter, colorlin);	
}

double Matrix::N_v_l_outside(const bool & colorlin,const int & Xiter, const int & Yiter)
{
	norm_assert(Xiter, Yiter, colorlin);
	this->double_fu = &Matrix::N_v_l_inside;
	return mode_switcher(Xiter, Yiter, colorlin);
}

double Matrix::N_v_e_outside(const bool & colorlin, const int & Xiter, const int & Yiter)
{
	norm_assert(Xiter, Yiter, colorlin);
	this->double_fu = &Matrix::N_v_e_inside;
	return mode_switcher(Xiter, Yiter, colorlin);
}

double Matrix::N_v_m_inside(const int & Xiter, const int & Yiter, int & size)
{
	double max = 0;
	for (int i = 0; i < ((size == 1) ? sizeX : sizeY); i++) {
		if(max < abs(mtrx[Yiter*sizeX + Xiter + size * i]))
		max = abs(mtrx[Yiter*sizeX + Xiter + size * i]);
	}
	return max; 
}

double Matrix::N_v_l_inside(const int & Xiter, const int & Yiter, int & size)
{
	int l = 0;
	for (int i = 0; i < ((size == 1) ? sizeX : sizeY); i++) {
		l+=abs(mtrx[Yiter*sizeX + Xiter + size * i]);
	}
	return l;
}

double Matrix::N_v_e_inside(const int & Xiter, const int & Yiter, int & size)
{
	int E = 0;
	for (int i = 0; i < ((size == 1 )? sizeX: sizeY); i++) {
		E += abs(mtrx[Yiter*sizeX + Xiter + size * i])*abs(mtrx[Yiter*sizeX + Xiter + size * i]);
	}
	double buf = pow(E, (1.0 / 2.0));
	return buf;
}

double Matrix::N_m_m()
{
	double max = 0;
	for (int i = 0; i < sizeX; i++) {
		double buf = 0;
		for (int j = 0; j < sizeY; j++)
			buf += abs(mtrx[j*sizeX + i]);
		max = max < buf ? buf : max;
	}
	return max;
}

double Matrix::N_m_c()
{
	double max = 0;
	for (int i = 0; i < sizeY; i++) {
		double buf = 0;
		for (int j = 0; j < sizeX; j++)
			buf += abs(mtrx[i*sizeX + j]);
		max = max < buf ? buf : max;
	}
	return max;
}
double Matrix::N_m_e()
{
	double sq = 0;
	for (int i = 0; i < sizeY; i++) {
		for (int j = 0; j < sizeX; j++) {
			sq += mtrx[i* sizeX + j] * mtrx[i* sizeX + j];
		}
	}
	return sqrt(sq);
}

double Matrix::N_m_s()
{
	double E = 0.001;
	double sq = 0;

	vector <double> x(sizeX, 0), y(sizeX, 0), z(sizeX, 0);
	x[0] = 1.0;
	Matrix A = T() * (*this);
	Matrix X{ x, sizeX, 1 };
	Matrix Y{ y,  sizeX, 1 };
	Matrix Z{ z,  sizeX, 1, };
	for (;;) {
		Y =  A * X;
		sq = Y.N_v_m_outside();
		Z = Y * ( 1.0 / sq);
		if ((Z.sum(X*(-1.0))).N_v_m_outside() < E)
			break;
		X = Z;
	}
    return sqrt(sq);
}

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
		Matrix b{ sizeX, sizeY };
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

Matrix  Matrix::operator*(const double & skalyar)
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
	Matrix a{ poryadok, poryadok };
	a.mtrx.resize(poryadok* poryadok);
	for (int i = 0; i < poryadok; i++) {
		for (int j = 0; j < poryadok; j++) {
			if (i == j)a.mtrx[i*poryadok + j] = 1;
			else
				a.mtrx[i*poryadok + j] = 0;

		}
	}
	return a;
}

double Matrix::tr()
{
	int i = 0, j = 0;

	if (sizeY == sizeX) {
		double  to_return = 0;
		for (; i < sizeY; i++) {
			j = i;
			to_return += mtrx[i*sizeX + j];

		}
		return to_return;
	}
	else return 0;
}

istream & operator>>(istream & in, Matrix & m)
{
	cout << "Enter data of the matrix: " << endl;
	for (int i = 0; i < m.sizeY; i++)
		for (int j = 0; j < m.sizeX; j++)
		{
			cout << "Enter i:" << i << " j: " << j;
			int buf = 0;
			m.Enter_check(buf);
			m.mtrx[i * m.sizeX + j] = buf;
		}
	return in;
}

ostream & operator<<(ostream & out, const  Matrix & m)
{
	cout << "Matrix: " << endl;
	for (int i = 0; i < m.sizeY; i++) {
		for (int j = 0; j < m.sizeX; j++)
		{
			cout << m.mtrx[i * m.sizeX + j] << "\t";
		}
		cout << endl;
	}

	return out;
}


