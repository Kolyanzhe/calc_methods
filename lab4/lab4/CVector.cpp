#include <vector>
#include <iostream>
#include "CVector.h"
#include <cmath>
#include "Header.h"
using namespace std;



vector<typenam> operator+(const vector<typenam>& a1, const vector<typenam>& a2) {
	vector<typenam> a3(a1);
	for (int i = 0; i < a3.size(); ++i) {
		a3[i] += a2[i];
	}
	return a3;
}

vector<typenam> operator-(const vector<typenam>& a1, const vector<typenam>& a2) {
	vector<typenam> a3(a1);
	for (int i = 0; i < a3.size(); ++i) {
		a3[i] -= a2[i];
	}
	return a3;
}

vector<typenam> operator*(const vector<typenam>& a, typenam x) {
	vector<typenam> result(a);
	for (int i = 0; i < a.size(); ++i)
		result[i] *= x;
	return result;
}

vector<typenam> operator*(typenam x, const vector<typenam>& a) {
	vector<typenam> result(a);
	for (int i = 0; i < a.size(); ++i)
		result[i] *= x;
	return result;
}

vector<typenam> operator/(const vector<typenam>& a, typenam x) {
	vector<typenam> result(a);
	for (int i = 0; i < a.size(); ++i)
		result[i] /= x;
	return result;
}

typenam operator*(const vector<typenam>& a1, const vector<typenam>& a2) {
	typenam result = 0.;
	for (int i = 0; i < a1.size(); ++i)
		result += a1[i] * a2[i];
	return result;
}

vector<typenam> operator+=(vector<typenam>& a1, const vector<typenam>& a2) {
	for (int i = 0; i < a1.size(); ++i)
		a1[i] += a2[i];
	return a1;
}

vector<typenam> operator-=(vector<typenam>& a1, const vector<typenam>& a2) {
	for (int i = 0; i < a1.size(); ++i)
		a1[i] += a2[i];
	return a1;
}

vector<typenam> operator*=(vector<typenam>& a, typenam x) {
	for (int i = 0; i < a.size(); ++i)
		a[i] *= x;
	return a;
}

vector<typenam> operator/=(vector<typenam>& a, typenam x) {
	for (int i = 0; i < a.size(); ++i)
		a[i] /= x;
	return a;
}
vector<typenam> normalize(const vector<typenam>& a) {
	vector<typenam> an(a);
	typenam norm = norm2(a);
	for (int i = 0; i < an.size(); ++i)
		an[i] /= norm;
	return an;
}

typenam sqr(typenam a) {
	return a * a;
}

typenam cube(typenam a) {
	return a * a * a;
}