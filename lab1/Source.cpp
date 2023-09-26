#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "Header.h"
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

vector <typenam> product(const vector<vector<typenam>> data, const vector<typenam>& b, int n) { // Óìíîæåíèå ìàòðèöû íà ñòîëáåö
    vector <typenam> result(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += data[i][j] * b[j];
        }
    }
    return result;
}

vector<typenam> qrDecomposition(vector< vector<typenam>>& A, vector< vector<typenam>>& R, vector< vector<typenam>>& Q, vector<typenam> b) {
    size_t n = A.size();
    R = A; // Копируем исходную матрицу A в R

    // Инициализируем матрицу Q как единичную матрицу
    Q = vector< vector<typenam>>(n, vector<typenam>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        Q[i][i] = 1.0;
    }

    for (size_t j = 0; j < n; ++j) {
        for (size_t i = n - 1; i >= j + 1; --i) { // Изменено условие
            typenam a = R[i - 1][j];
            typenam b = R[i][j];
            typenam r = sqrt(a * a + b * b);

            // Вычисляем c и s
            typenam c = a / r;
            typenam s = -b / r;

            // Обновляем матрицу R
            for (size_t k = 0; k < n; ++k) {
                typenam temp = c * R[i - 1][k] - s * R[i][k];
                R[i][k] = s * R[i - 1][k] + c * R[i][k];
                R[i - 1][k] = temp;

            }

            // Обновляем матрицу Q
            for (size_t k = 0; k < n; ++k) {
                typenam temp = c * Q[i - 1][k] - s * Q[i][k];
                Q[i][k] = s * Q[i - 1][k] + c * Q[i][k];
                Q[i - 1][k] = temp;
            }
        }
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            R[i][j] = 0.0;
        }
    }
    vector<vector<typenam>>Q_T(n, vector<typenam>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            Q_T[i][j] = Q[j][i];
        }
    }
    vector<typenam> b1(n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            b1[i] += Q_T[j][i] * b[j];
        }
    }

    vector <typenam> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b1[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= R[i][j] * x[j];
        }
        x[i] /= R[i][i];
    }
        return x;

}



