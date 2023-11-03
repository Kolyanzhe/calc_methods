#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "Header.h"
#include "Matrix.h"
#include "CVector.h"
using namespace std;


typenam norm1(const vector<typenam>& a) { //Îêòàýäðè÷åñêàÿ íîðìà 
    typenam result = 0;
    for (int i = 0; i < a.size(); ++i)
        result += fabs(a[i]);
    return result;
}

typenam norm2(const vector<typenam>& a) { // Åâêëèäîâà íîðìà 
    typenam result = 0;
    for (int i = 0; i < a.size(); ++i)
        result += (a[i] * a[i]);
    return sqrt(result);
}

typenam normInf(const vector<typenam>& a) { // Êóáè÷åñêàÿ íîðìà 
    typenam max = 0;
    for (int i = 0; i < a.size(); ++i)
        if (fabs(a[i]) > max)
            max = fabs(a[i]);
    return max;
}

vector <typenam> productmv(const vector<vector<typenam>> data, const vector<typenam>& b, int n) {
    // Óìíîæåíèå ìàòðèöû íà ñòîëáåö
    vector <typenam> result(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += data[i][j] * b[j];
        }
    }
    return result;
}

// Óìíîæåíèå ìàòðèö
vector<vector<typenam>> product(const vector<vector<typenam>>& a, const vector<vector<typenam>>& b, int n) {
    vector<typenam> r(n, 0);
    vector<vector<typenam>> cdot(n, r);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                cdot[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return cdot;
}

vector<vector<typenam>> inversematrix(const vector<vector<typenam>>& data, int n) {
    vector <typenam> r(n, 0);
    vector<vector<typenam>> reversematrix(n, r);
    Matrix A(n, data);
    vector<vector<typenam>> Q = A.qrDec().second;
    vector<vector<typenam>> R = A.qrDec().first;
    vector<typenam> v(n, 0);
    vector<typenam> e(n, 0);
    for (int i = 0; i < n; ++i) {
        e[i] = 1;
        v = A.qrDecomposition(Q, R, e);
        for (int j = 0; j < n; j++) {
            reversematrix[j][i] = v[j];
        }
        v = r;
        e = r;
    }
    return reversematrix;
}
