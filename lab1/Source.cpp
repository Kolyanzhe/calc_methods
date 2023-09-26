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

void print(vector<vector<typenam>>& matrix) {
    for (const auto& row : matrix) {
        for (const typenam& element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}


