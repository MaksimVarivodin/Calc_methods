#include <iostream>
#include <cmath>
using namespace std;
// функция формирования массива дельт
double*& delta_mas(double *& init_mas, int &size) {
	
	double* del = new double[size];
	double* buf = new double[size - 1];
	int SIZE = size;
	del[0] = init_mas[0];// getting delta0 y0
	cout << endl;
	cout << "y0: "<< endl;
	for (int i = SIZE - 1; i >= 0; i--)
		cout << init_mas[i]<< " ";
	cout<< endl << "delta massivs:" << endl;
	for (int i = 1; i < SIZE; i++) {
		cout << "del " << i << " y0" << endl;
		for (int j = size - 1; j > 0; j--) {
			buf[j - 1] = init_mas[j] - init_mas[j - 1];
			cout << buf[j - 1]<< " ";
		} 
		cout << endl;
		delete[] init_mas;
		del[i] = buf[0];// to-return delta1-4 y0
		init_mas = buf;
		size--;
		buf = new double[size - 1];
	}
	delete[] init_mas;
	return del;
}
double* pol(double *masX, double* del, const int &size = 5) {
	double *res = new double[size];
	double *masXP = new double[size];
	// задаем массив иксов точек
	for (int i = 0; i < size; i++)
		masXP[i] = masX[i]+ 0.6;
	// считаем шаг
	double h = masX[1] - masX[0];
	cout << endl << "y: "<< endl;
	// общий цикл для игриков
	for (int k = 0; k < size; k++) {
		res[k] = 0;
		// расчет одного игрика
		for (int i = 0; i < size; i++) {
			double mult = del[i];// дельта
			for (int j = 0; j < i; j++) {
				mult *= (masXP[k] - masX[j]);
				mult /= (j+ 1)*h;
			}
								
			res[k] += mult;
		}
	}
	for (int k = 0; k < size; k++) {
		
		cout << res[k]<< " " << endl;
	}		
	return res;
}
double *& NewT() {	
	int size = 5, size0 = size;
	double *masX = new double [size] {0, 1.2, 2.4, 3.6, 4.8};
	cout << "    0  1  2  3  4" << endl<< "x: ";
	for (int i = 0; i < size; i++)
		cout << masX[i] << " ";

	double *masY = new double[size] {8.6, -12.2, -2.8, 2.7, -5.3};
	double* del = delta_mas(masY, size0);
	cout  << endl << "d: ";
	for (int i = 0; i < size; i++)
		cout << del[i] << " ";
	double *masYP = pol(masX, del);	

	delete[] masX;	
	return masYP;
}
int main()
{
	double* masYp = NewT();
	system("pause");
	delete[] masYp;
}