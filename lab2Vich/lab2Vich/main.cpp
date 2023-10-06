#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "Header.h"
#include "Matrix.h"
#include "CVector.h"
using namespace std; 

int main() {
    vector<typenam> x1 = { 5,-7,12,4 };
    vector<typenam> x2 = { 10,-10,12,4 };
    ifstream file;   //файл с размерностью, левой частью и столбцом правой части 
    file.open("data1.txt");
    int variant = 4;
    int n;  //размерность матрицы 
    file >> n; // Это размерность 
    typenam s;
    vector<typenam> v(n, 0);
    vector<vector<typenam>> data(n, v);
    for (int i = 0; i < n; i++) {  //считываем матрицу
        for (int j = 0; j < n; j++) {
            file >> s;
            data[i][j] = s;
        }
    }
    Matrix A(n, data);

    vector<typenam> b(n, 0);
    for (int i = 0; i < n; i++) {  //считываем столбец правой части   
        file >> s;
        b[i] = s;
    }

    cout << "Matrix A:" << endl; // Вывод матрицы
    A.print();
    cout << "Vector b:" << endl;  // Вывод столбца (заданного)
    for (int i = 0; i < n; i++) {
        cout << "b[" << i << "] = " << b[i] << endl;
    }
    
    cout << "*** JACOBI METHOD ***" << endl << endl;
    vector<typenam> x = A.Jacobi(b);
    cout << "Solution: " << endl;
    for (size_t i = 0; i < n; ++i) {
        cout << "x[" << i + 1 << "] = " << x[i] << endl;
    }
    cout << "Norma nevyazki: " << norm1(productmv(data,x,n)-b) << endl;
    cout << "error: " << norm2(x2-x) << endl;

    cout << endl << "*** ZEIDEL METHOD ***" << endl << endl;
    x = A.Zeidel(b);
    cout << "Solution: " << endl;
    for (size_t i = 0; i < n; ++i) {
        cout << "x[" << i + 1 << "] = " << x[i] << endl;
    }
    cout << "Norma nevyazki: " << norm1(productmv(data, x, n) - b) << endl;
    cout << "error: " << norm2(x2 - x) << endl;

    //Большая система Зейдель метод
    cout << endl << "SLAE n = 204, Ziedel Method" << endl << endl;
    int N = 200 + variant;
    vector<typenam> a_slae(N, 1);
    a_slae[0] = 0;
    vector<typenam> b_slae(N, 4);
    vector<typenam> c_slae(N, 1);
    c_slae[N - 1] = 0;
    vector<typenam> d_slae(N, 0);
    for (size_t i = 0; i < N; ++i) { d_slae[i] = 10 - 2 * (i % 2); }
    d_slae[0] = 6;
    d_slae[N - 1] = 9.0 - 3.0 * ((N) % 2);
    vector<typenam> X = BigSlaeZeidel(a_slae, b_slae, c_slae, d_slae, N);
    vector<typenam> X_sol(N, 0);
    for (size_t i = 0; i < N; ++i) {  //Известное решение 
        X_sol[i] = 2 - (i % 2);
    }
    for (size_t i = 0; i < N; ++i) {
        cout << "x[" << i + 1 << "] = " << X[i] << endl;
    }
    cout << norm2(X - X_sol) << endl;



}








