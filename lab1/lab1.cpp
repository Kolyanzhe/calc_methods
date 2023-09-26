#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "Header.h"
#include "CVector.h"
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
    //Matrix A(n, n, data);

    vector<typenam> b(n, 0);
    for (int i = 0; i < n; i++) {  //ñ÷èòûâàåì ñòîëáåö ïðàâîé ÷àñòè   
        file >> s;
        b[i] = s;
    }

    //cout << "Matrix A:" << endl; // Âûâîä ìàòðèöû
    //A.print();
    //cout << "Vector b:" << endl;  // Âûâîä ñòîëáöà (çàäàííîãî)
    //for (int i = 0; i < n; i++) {
    //    cout << "b[" << i << "] = " << b[i] << endl;
    //}
    //vector<typenam> x = A.solveGaussPartialPivoting(b);
    //cout << "Solution:" << endl;
    //for (int i = 0; i < 4; ++i) {
    //    cout << "x[" << i + 1 << "] = " << x[i] << endl;
    //}
    //cout << "Vector b1=Ax: " << endl;
    //vector <typenam> b1 = product(data, x, n);
    //for (int i = 0; i < n; i++) {      // Âûâîä âåêòîðà b1=A*x
    //    cout << "b1[" << i + 1 << "] = " << b1[i] << endl;
    //}
    //vector <typenam> b_b1(b - b1);  // Âåêòîð íåâÿçêè 
    //cout << "b-b1: " << endl;
    //for (int i = 0; i < n; i++) {      // Âûâîä âåêòîðà íåâÿçêè
    //    cout << "(b-b1)[" << i + 1 << "] = " << b_b1[i] << endl;
    //}
    //cout << "Norma nevyazki: " << norm1(b_b1) << endl; 

    // Пример исходной матрицы A
    vector<vector<double>> A = data;
    
    n = A.size();

    // Создаем матрицы Q и R
    vector<vector<double>> Q, R;
    vector<typenam> x;
    // Выполняем QR-разложение
    x = qrDecomposition(A, R, Q, b);
    for (size_t i = 0; i < x.size(); i++) {
        cout << x[i] << endl;
    }
    return 0;
}


