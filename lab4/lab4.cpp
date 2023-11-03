#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "Header.h"

using namespace std;

int main() {
    ifstream file;   //ôàéë ñ ðàçìåðíîñòüþ, ëåâîé ÷àñòüþ è ñòîëáöîì ïðàâîé ÷àñòè 
    file.open("data1.txt");

    int n;  //ðàçìåðíîñòü ìàòðèöû 
    file >> n; // Ýòî ðàçìåðíîñòü 
    typenam s;
    vector<typenam> v(n, 0);
    vector<vector<typenam>> data(n, v);
    for (int i = 0; i < n; i++) {  //ñ÷èòûâàåì ìàòðèöó
        for (int j = 0; j < n; j++) {
            file >> s;
            data[i][j] = s;
        }
    }
    Matrix A(n, data);

    A.print();
    A.QR_lambda_sdvig();
}

