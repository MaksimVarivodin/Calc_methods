#include "Matrix4_cl.h"


Matrix Matrix::X(Matrix d)
{
	 double max;
	 Matrix x{sizeX, 1};
	 Matrix a{ *this };
	 Matrix b{d};
	int k, index;
	const double eps = 0.00001;  // точность
	
	k = 0;
	while (k < sizeX)
	{
		// Поиск строки с максимальным a[i][k]
		max = abs(a.mtrx[k*sizeX + k]);
		index = k;
		for (int i = k + 1; i < sizeX; i++)
		{
			if (abs(a.mtrx[i*sizeX + k]) > max)
			{
				max = abs(a.mtrx[i*sizeX + k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			cout << "Решение получить невозможно из-за нулевого столбца ";
			cout << index << " матрицы A" << endl;
			return Matrix{0, 0};
		}
		for (int j = 0; j < sizeX; j++)
		{
			double temp = a.mtrx[k*sizeX + j];
			a.mtrx[k*sizeX + j] = a.mtrx[index*sizeX + j];
			a.mtrx[index*sizeX + j] = temp;
		}
		double temp = b.mtrx[k];
		b.mtrx[k] = b.mtrx[index];
		b.mtrx[index] = temp;
		// Нормализация уравнений
		for (int i = k; i < sizeX; i++)
		{
			double temp = a.mtrx[i*sizeX + k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < sizeX; j++)
				a.mtrx[i*sizeX + j] = a.mtrx[i*sizeX + j] / temp;
			b.mtrx[i] = b.mtrx[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < sizeX; j++)
				a.mtrx[i*sizeX + j] = a.mtrx[i*sizeX + j] - a.mtrx[k*sizeX + j];
			b.mtrx[i] = b.mtrx[i] - b.mtrx[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = sizeX - 1; k >= 0; k--)
	{
		x.mtrx[k] = b.mtrx[k];
		for (int i = 0; i < k; i++)
			b.mtrx[i] = b.mtrx[i] - a.mtrx[i* sizeX + k] * b.mtrx[k];
	}
	return x;
	
}
// Условие окончания
bool Matrix::converge(Matrix xk, Matrix xkp, double eps)
{
	double norm = 0;
	for (int i = 0; i < xk.sizeX; i++)
		norm += (xk.mtrx[i] - xkp.mtrx[i]) * (xk.mtrx[i] - xkp.mtrx[i]);
	return (sqrt(norm) < eps);
}

double Matrix::okr(double x, double eps)
{
	int i = 0;
	double neweps = eps;
	while (neweps < 1)
	{
		i++;
		neweps *= 10;
	}
	int okr = pow(double(10), i);
	x = int(x * okr + 0.5) / double(okr);

	return x;
}

bool Matrix:: diagonal(Matrix a)
{
	int i, j, k = 1;
	double sum;
	for (i = 0; i < a.sizeX; i++) {
		sum = 0;
		for (j = 0; j < a.sizeX; j++) sum += abs(a.mtrx[i*a.sizeX + j]);
		sum -= abs(a.mtrx[i*a.sizeX + i]);
		if (sum > a.mtrx[i*a.sizeX + i])
		{
			k = 0;
		}
	}

	return (k == 1);

}
Matrix Matrix::Zejdel(Matrix b)
{
	Matrix x{ sizeX, 1 };
	Matrix a{* this};
	Matrix p{sizeX, 1};
	int m = 0;
	double eps = 0.0001;
	if (diagonal(a)) {
		do
		{
			for (int i = 0; i < sizeX; i++)
				p.mtrx[i] = x.mtrx[i];
			for (int i = 0; i < sizeX; i++)
			{
				double var = 0;
				for (int j = 0; j < sizeX; j++)
					if (j != i) var += (a.mtrx[i* sizeX + j] * x.mtrx[j]);

				x.mtrx[i] = (b.mtrx[i] - var) / a.mtrx[i* sizeX + i];
			}
			m++;
		} while (!converge(x, p, eps));




		for (int i = 0; i < sizeX; i++) x.mtrx[i] = okr(x.mtrx[i], eps);
		cout << "Iterations: " << m << endl;
	}
	else {
		cout << "Не выполняется преобладание диагоналей" << endl;
	}
	return x;
}
double Matrix::cond1()
{
	return N_m_m() * (m_minus_1().N_m_m());
}

double Matrix::cond2()
{
	return N_m_s() * (m_minus_1().N_m_s());
}

double Matrix::condI()
{
	return N_m_c() * (m_minus_1().N_m_c());
}

double Matrix::condE()
{
	return N_m_e() * (m_minus_1().N_m_e());
}
Matrix Matrix::m_minus_1()
{
	Matrix buf{ *this };
	return buf * (1.0 / det_calc());
}



Matrix Matrix::U_calc()
{
	vector<double> E; // заготовка для матриц М (единичная матрица)
	E.resize(sizeX * sizeY);
	for (int i = 0; i < sizeY; i++)
		for (int j = 0; j < sizeX; j++)
			E[i*sizeX + j] = (i == j ? 1 : 0);
	// для первой итерации нужно чтобы М1* А
	Matrix buf{ *this };

	for (int i = 0; i < sizeX - 1; i++) {
		Matrix Mn{ E, sizeX, sizeY };
		for (int j = i + 1; j < sizeY; j++) {
			if (buf.mtrx[i*sizeX + i] == 0) {
				buf.line_swap(i, i);
			}
			Mn.mtrx[j*sizeX + i] = (buf.mtrx[j*sizeX + i] / buf.mtrx[i*sizeX + i]) * (-1);
		}
		buf = Mn * buf;
	}
	return buf;
}

Matrix Matrix::L_calc()
{
	vector<double> E; // заготовка для матриц М (единичная матрица)
	E.resize(sizeX * sizeY);
	for (int i = 0; i < sizeY; i++)
		for (int j = 0; j < sizeX; j++)
			E[i*sizeX + j] = (i == j ? 1 : 0);
	vector <Matrix> mas;
	vector<Matrix>::iterator it = mas.begin();
	mas.insert(it, sizeX - 1, Matrix{ E, sizeX, sizeY });
	// для первой итерации нужно чтобы М1* А
	Matrix buf{ *this };

	for (int i = 0; i < sizeX - 1; i++) {

		for (int j = i + 1; j < sizeY; j++) {
			if (buf.mtrx[i*sizeX + i] == 0) {
				buf.line_swap(i, i);
			}
			mas[i].mtrx[j*sizeX + i] = (buf.mtrx[j*sizeX + i] / buf.mtrx[i*sizeX + i]) * (-1);
		}

		buf = mas[i] * buf;
	}
	Matrix L{
		E, sizeX, sizeY
	};
	for (int i = 0; i < sizeX - 1; i++) {
		for (int j = i + 1; j < sizeY; j++) {
			L.mtrx[j*sizeX + i] = mas[i].mtrx[j*sizeX + i] * (-1);
		}
	}
	return L;
}

double Matrix::det_calc()
{
	Matrix buf = U_calc();
	double n = 1;
	for (int i = 0; i < sizeX; i++)
		n *= buf.mtrx[i*sizeX + i];
	return n;
}
void Matrix::line_swap(const int & i, const int & j)
{
	for (int j0 = j + 1; j0 < sizeY; j0++)
	{
		if (mtrx[j0* sizeX + i] != 0) {
			Matrix buf{ sizeX, 1 };
			for (int i0 = 0; i0 < sizeX; i0++) {
				buf.mtrx[i0] = mtrx[j0*sizeX + i0];
				mtrx[j0*sizeX + i0] = mtrx[j*sizeX + i0];
				mtrx[j*sizeX + i0] = buf.mtrx[i0];
			}
			return;
		}

	}

	for (int i0 = i + 1; i0 < sizeX; i0++)
	{
		if (mtrx[j* sizeX + i0] != 0) {
			Matrix buf{ 1, sizeY };
			for (int j0 = 0; j0 < sizeY; j0++) {
				buf.mtrx[j0] = mtrx[j0 * sizeX + i0];
				mtrx[j0*sizeX + i0] = mtrx[j0*sizeX + i];
				mtrx[j0*sizeX + i] = buf.mtrx[j0];
			}
			return;
		}

	}

}

void Matrix::norm_assert(const int & Xiter, const int & Yiter, const bool & colorlin)
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

double Matrix::N_v_l_outside(const bool & colorlin, const int & Xiter, const int & Yiter)
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
		if (max < abs(mtrx[Yiter*sizeX + Xiter + size * i]))
			max = abs(mtrx[Yiter*sizeX + Xiter + size * i]);
	}
	return max;
}

double Matrix::N_v_l_inside(const int & Xiter, const int & Yiter, int & size)
{
	int l = 0;
	for (int i = 0; i < ((size == 1) ? sizeX : sizeY); i++) {
		l += abs(mtrx[Yiter*sizeX + Xiter + size * i]);
	}
	return l;
}

double Matrix::N_v_e_inside(const int & Xiter, const int & Yiter, int & size)
{
	int E = 0;
	for (int i = 0; i < ((size == 1) ? sizeX : sizeY); i++) {
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
		Y = A * X;
		sq = Y.N_v_m_outside();
		Z = Y * (1.0 / sq);
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





