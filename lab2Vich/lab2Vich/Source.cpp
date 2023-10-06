#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "Header.h"
#include "Matrix.h"
#include "CVector.h"
using namespace std;


typenam norm1(const vector<typenam>& a) { //Октаэдрическая норма 
    typenam result = 0;
    for (int i = 0; i < a.size(); ++i)
        result += fabs(a[i]);
    return result;
}

typenam norm2(const vector<typenam>& a) { // Евклидова норма 
    typenam result = 0;
    for (int i = 0; i < a.size(); ++i)
        result += (a[i] * a[i]);
    return sqrt(result);
}

typenam normInf(const vector<typenam>& a) { // Кубическая норма 
    typenam max = 0;
    for (int i = 0; i < a.size(); ++i)
        if (fabs(a[i]) > max)
            max = fabs(a[i]);
    return max;
}

vector <typenam> productmv(const vector<vector<typenam>> data, const vector<typenam>& b, int n) { 
    // Умножение матрицы на столбец
    vector <typenam> result(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += data[i][j] * b[j];
        }
    }
    return result;
}

// Умножение матриц
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

// Здесь обратную матрицу ищем
vector<vector<typenam>> inversematrix(const vector<vector<typenam>>& data, int n) {
    vector <typenam> r(n, 0);
    vector<vector<typenam>> reversematrix(n, r);
    Matrix A( n, data);
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


vector<typenam> BigSlaeZeidel(const vector<typenam>& a, const vector<typenam>& b, const vector<typenam>& c, const vector<typenam>& d, int N) {
    vector<typenam> nul(N, 0);
    vector<typenam> x_0(N, 1);
    vector<typenam>x_k(N, 0);
    vector<typenam> x_k_x_0 = x_0;
    typenam Usum = 0;
    typenam Lsum = 0;
    //Для остановы 
    vector<typenam> AX(N, 0);
    typenam nevyazka = 0;

    size_t iter_culc = 1;
    while ((norm1(x_k_x_0) / (norm1(x_k) + 0.0000000001) >= epsilon) || (norm2(AX - d) >= epsilon)) {
        AX = nul;
        for (size_t i = 1; i < N - 1; ++i) {

            Usum += a[i] * x_0[i - 1];
            Lsum += c[i] * x_k[i + 1];
            x_k[i] = (d[i] - Usum - Lsum) / b[i];

            Usum = 0;
            Lsum = 0;
        }
        x_k[0] = (d[0] - c[0] * x_0[0]) / b[0];
        x_k[N - 1] = (d[N - 1] - a[N - 1]*x_0[N - 1]) / b[N - 1];
        for (size_t k = 0; k < N; ++k) {
            x_k_x_0[k] = x_0[k] - x_k[k];
        }
        x_0 = x_k;
        ++iter_culc;
        if (iter_culc > 100000) {
            cout << "TOO MANY ITERATIONS" << endl;
            return x_k;
        }
        for (size_t i = 1; i < N-1; ++i) {
            AX[i] += a[i] * x_0[i - 1] + b[i] * x_0[i] + c[i] * x_0[i + 1];
        }
        AX[0] = b[1] * x_0[1] + c[1] * x_0[2];
        AX[N - 1] = a[N - 1] * x_0[N - 2] + b[N - 1] * x_0[N - 1];
    }
    return x_k;
}













