#include <NTL/ZZ.h>
#include <iostream>
#include <ctime>
#include <cassert>

using namespace NTL;
using namespace std;

ZZ my_jacobi(ZZ a, ZZ b)
{
	if (GCD(a, b) != 1)
		return conv<ZZ>(0);
	ZZ r(1);
	if (a < 0)
	{
		a = -a;
		if (b % 4 == 3)
			r = -r;
	}
step4:
	ZZ t(0);
	while (a % 2 == 0)
	{
		++t;
		a /= 2;
	}
	if (t % 2 == 1)
		if (b % 8 == 3 || b % 8 == 5)
			r = -r;

	if (a % 4 == 3 && b % 4 == 3)
		r = -r;
	ZZ c = a;
	a = b % c;
	b = c;
	if (a != 0)
		goto step4;
	else
		return r;
}

ZZ LCPRNG(const ZZ& n)
{
	assert(n > 0);
	ZZ a = (ZZ)rand() % n + 2; // (a >= 2)
	ZZ x0 = (ZZ)rand() % n; // (x0 >= 0)
	ZZ c = (ZZ)rand() % n; // (c >= 0)

	for (int i = 0; i < 5; i++)
	{
		a = a * a;
		x0 = x0 * x0;
		c = c * c;
	}

	a = a < n ? a : a % n; // (a < n)
	x0 = x0 < n ? x0 : x0 % n; // (x0 < n)
	c = c < n ? c : c % n; // (c < n)

	return (a * x0 + c) % n;
}

ZZ ModExp(const ZZ& b, const ZZ& e, const ZZ& m) //b = base, e = exponent, m = modular
{
	long k = NumBits(e);

	ZZ res;
	res = 1;

	for (long i = k - 1; i >= 0; i--) {
		res = (res*res) % m;
		if (bit(e, i) == 1) res = (res*b) % m;
	}

	return res;
}

long solovay_strassen(const ZZ& n, const long& rounds = 10)
{
	srand(time(NULL));
	for (int i = 0; i < rounds; ++i)
	{
		long size = 2 + rand() % n.size(); // choose a random size of a
		ZZ a = LCPRNG(n);
		auto x = my_jacobi(a, n);
		if (x == 0)
			return 0;
		ZZ y = ModExp(a, (n - 1) / 2, n);
		if ((x - y) % n != 0)
			return 0;
	}
	return 1;
}

ZZ gen_prime(const long & size) // generate a random prime number
{
	ZZ n = RandomLen_ZZ(size);
	if (n % 2 == 0)
		++n;
	while (!solovay_strassen(n))
		n += 2;
	return n;
}


int main()
{
	for (int i = 0; i < 30; ++i)
	{
		ZZ n = RandomLen_ZZ(128); // generate arbitrary number
		cout << solovay_strassen(n) << ' '; // most likely 0
	}
	cout << endl;

	for (int i = 0; i < 30; ++i)
	{
		ZZ n = RandomPrime_ZZ(128); // generate a prime number
		cout << solovay_strassen(n) << ' '; // should be 1
	}
	cout << endl;

	ZZ n = conv<ZZ>("3490529510847650949147849619903898133417764638493387843990820577"); // it's a prime number
	cout << solovay_strassen(n) << endl; // must be 1

	ZZ n1 = conv<ZZ>("112505380450296606970338459629988782604252033209350010888227147338120001"); // it is not a prime number (Carmichael number)
	cout << solovay_strassen(n1) << endl; // must be 0
}
