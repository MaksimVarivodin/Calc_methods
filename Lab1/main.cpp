#include "lagr_pol_cl.h"
int  main(int argc) {
	setlocale(0, "");
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	NEWT A;
	cin >> A;
	A.polynom();
	cout << A;
	A.~NEWT();
	system("pause");
}