#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

#define typenam double
const typenam epsilon = 0.0001;

#include "Matrix.h"
#include"CVector.h"
typenam norm1(const vector<typenam>& a); //Октаэдрическая норма
typenam norm2(const vector<typenam>& a); // Евклидова норма 
typenam normInf(const vector<typenam>& a); // Кубическая норма
vector <typenam> productmv(const vector<vector<typenam>> data, const vector<typenam>& b, int n);  // Умножение матрицы на столбец
vector<vector<typenam>> product(const vector<vector<typenam>>& a, const vector<vector<typenam>>& b, int n); //Умножение матриц
vector<vector<typenam>> inversematrix(const vector<vector<typenam>>& data, int n);
vector<typenam> BigSlaeZeidel(const vector<typenam>& a, const vector<typenam>& b, const vector<typenam>& c, const vector<typenam>& d, int N);


class Matrix {
private:
    int n;
    vector<vector<typenam>> data;

public:
    Matrix(int n, vector<vector<typenam>> data) {
        this->n = n;
        this->data = data;
    }



    void print() const {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                cout << data[i][j] << "    ";
            }
            cout << endl;
        }
    }
    int size() const { return n; }

    vector<vector<typenam>> get_data() const { return data; }

    Matrix cdot(const Matrix a) {
        vector<vector<typenam>> cdot = data;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    cdot[i][j] += cdot[i][k] * a.get_data()[k][j];
                }
            }
        }
        Matrix CDOT(n, data);
        return CDOT;
    }


    //Матрица + Матрица
    Matrix operator+=(const Matrix& matr) {

        for (size_t row = 0; row < n; row++)
            for (size_t col = 0; col < n; col++)
                data[row][col] += matr.data[row][col];
        return *this;
    }

    //Матрица - Матрица
    Matrix operator-=(const Matrix& matr) {

        for (size_t row = 0; row < n; row++)
            for (size_t col = 0; col < n; col++)
                data[row][col] -= matr.data[row][col];
        return *this;
    }


    // Сложение матриц
    Matrix operator+(const Matrix& matr) const {
        Matrix result(*this);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                result.data[i][j] += matr.data[i][j];
        return result;
    }

    //Вычитание матриц
    Matrix operator-(const Matrix& matr) const {
        Matrix result(*this);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                result.data[i][j] -= matr.data[i][j];
        return result;
    }

    //Умножение матриц
    Matrix operator *(const Matrix& matr) const {
        vector<typenam> r(n, 0);
        vector<vector<typenam>> res(n, r);
        Matrix result(n, res);
        Matrix now(*this);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                for (size_t k = 0; k < n; ++k) {
                    result.data[i][j] += now.data[i][k] * matr.data[k][j];
                }
            }
        }
        return result;
    }

    //Умножение матрицы на число
    Matrix operator*(typenam k) const {
        Matrix result(*this);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                result.data[i][j] *= k;
        return result;
    }

    // Октаэдрическая матричная норма 
    typenam mnorm1() {
        typenam cond = 0;
        typenam max = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; j++) {
                max += fabs(data[j][i]);
            }
            if (cond < max) { cond = max; }
            max = 0;
        }
        return cond;
    }

    //Кубическая матричная норма
    typenam mnormInf() {
        typenam cond = 0;
        typenam max = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; j++) {
                max += fabs(data[i][j]);
            }
            if (cond < max) { cond = max; }
            max = 0;
        }
        return cond;
    }

    pair<vector< vector<typenam>>, vector<vector<typenam>>> qrDec() {
        vector<vector<typenam>> R(data);
        bool isUpperTriangular = true;
        for (size_t i = 1; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                if (R[i][j] != 0.0) {
                    isUpperTriangular = false;
                    break;
                }
            }
            if (!isUpperTriangular) {
                break;
            }
        }
        // Инициализируем матрицу Q как единичную матрицу
        vector< vector<typenam>>Q(n, vector<typenam>(n, 0.0));
        for (size_t i = 0; i < n; ++i) {
            Q[i][i] = 1.0;
        }
        if (!isUpperTriangular) {

            
            typenam a, c, s, buf;
            size_t ii, jj;
            for (size_t i = 0; i < n - 1; ++i) {

                for (size_t j = i + 1; j < n; ++j) {
                    ii = i;
                    jj = j;
                    if (fabs(R[i][j]) > 1e-9) {
                        a = sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
                        c = R[i][i] / a;
                        s = R[j][i] / a;
                        for (size_t k = i; k < n; ++k) {
                            buf = R[i][k];
                            R[i][k] = c * buf + s * R[j][k];
                            R[j][k] = -s * buf + c * R[j][k];
                        }
                        for (size_t k = 0; k < n; ++k) {
                            buf = Q[i][k];
                            Q[i][k] = c * buf + s * Q[j][k];
                            Q[j][k] = -s * buf + c * Q[j][k];

                        }

                        R[jj][ii] = 0;
                    }
                }
            }

        }
        return pair<vector< vector<typenam>>, vector<vector<typenam>>>(R, Q);
    }
    vector<typenam> qrDecomposition(vector<vector<typenam>>& QQ, vector<vector<typenam>>& RR, vector<typenam>& b) {      
        vector<vector<typenam>> R = RR;
        vector<vector<typenam>> Q = QQ;
        //Проверка верхнетреугольности
        bool isUpperTriangular = true;
        for (size_t i = 1; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                if (data[i][j] != 0.0) {
                    isUpperTriangular = false;
                    break;
                }
            }
            if (!isUpperTriangular) {
                break;
            }
        }
        vector <typenam> x(n);
        
        if (!isUpperTriangular) {
            //Транспонируем Q(можно и без этого)
            vector<vector<typenam>>Q_T(n, vector<typenam>(n, 0.0));
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    Q_T[i][j] = Q[j][i];
                }
            }
            vector<typenam> b1(n, 0.0);
            //Считаем b*
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    b1[i] += Q_T[j][i] * b[j];
                }
            }
            //Обратный ход
            for (int i = n - 1; i >= 0; --i) {
                if (fabs(R[i][i]) < 1e-9) {
                    cout << "Matrica virozhdena" << endl;
                    break;
                }
                x[i] = b1[i];
                for (int j = i + 1; j < n; ++j) {
                    x[i] -= R[i][j] * x[j];
                }
                x[i] /= R[i][i];
            }

            return x;
        }
        else {
            for (int i = n - 1; i >= 0; --i) {
                x[i] = b[i];
                for (int j = i + 1; j < n; ++j) {
                    x[i] -= R[i][j] * x[j];
                }
                x[i] /= R[i][i];
            }
            return x;
        }

    }



    //Метод Якоби

    vector<typenam> Jacobi(const vector<typenam>& b) { 
        vector<typenam> x_0 =  b; // Начальное приближение 
        vector<typenam> diag(n, 0); // Это диагональная матрица D=B (x=Bx+c)
        size_t iter_culc = 1;
        for (size_t i = 0; i < n; ++i) {
            if (fabs(data[i][i]) < 1e-9) {
                cout << "NULL ON DIAG" << endl;
                exit(1);
            }
            diag[i] = data[i][i];
        }
        vector<typenam> diag_inverse(n, 0);
        for (size_t i = 0; i < n; ++i) {
            diag_inverse[i] = 1/diag[i];
        }
        vector<vector<typenam>> Matrix_B = data;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i == j) { Matrix_B[i][i] = 0; }
                else Matrix_B[i][j] *= -diag_inverse[i];
            }
        }
        Matrix B(n, Matrix_B);
        typenam matr_norm_1 = B.mnorm1();
        typenam matr_norm_Inf = B.mnormInf();
        cout << "Matrix B norm1: " << matr_norm_1 << endl;
        cout << "Matrix B norm_inf: " << matr_norm_Inf << endl;
        vector<typenam>x_k(n, 0);
        vector<typenam> x_k_x_0 = x_0;
        typenam sum = 0;
        if (matr_norm_Inf >= 1) {
            cout << "||C||>1" << endl;
            exit(1);
        }
        while (norm1(x_k_x_0) >= (1- matr_norm_Inf)/ matr_norm_Inf *epsilon) {
           // x_k = x_0;
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    if (i != j) { sum += data[i][j] * x_0[j]; }
                }
                x_k[i] = (b[i] - sum) / data[i][i];
                sum = 0;
            }
            for (size_t k = 0; k < n; ++k) {
                x_k_x_0[k] = x_0[k] - x_k[k];
            }
            x_0 = x_k;
            ++iter_culc;
            if (iter_culc>1000) {
                cout << "TOO MANY ITERATIONS" << endl;
                exit(1);
            }
        }
       cout << "A priori estimate iterations: " <<int( log(epsilon) / log( matr_norm_Inf)) << endl;
       cout << "A posteriori estimate iterations: " << int(log(epsilon) / log(matr_norm_Inf / (1 - matr_norm_Inf))) << endl;
       cout << "Iterations: " << iter_culc << endl;
        return x_k;
    }

    //Метод Зейделя

    vector<typenam> Zeidel(const vector<typenam>& b) {
        vector<typenam> x_0 = b; // Начальное приближение 
        size_t iter_culc = 1;
        for (size_t i = 0; i < n; ++i) {
            if (fabs(data[i][i]) < 1e-9) {
                cout << "NULL ON DIAG" << endl;
                exit(1);
            }
        }
        vector<vector<typenam>> B = data;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i ; j < n; ++j) {
                B[i][j] = 0;
            }
        }
        Matrix A1(n, B);
        for (size_t i = 0; i < n; i++) {
            B[i][i] = data[i][i];// (A1+D)x^{k+1}+A2 x^k = f; A1+D = B
        }
        vector<vector<typenam>> A2 = data;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                A2[j][i] = 0;
            }
        }
        
        Matrix BB(n, B);
        Matrix AA2(n, A2);  // Это то же самое, что и A2, но Matrix! Нужно, чтобы 


        Matrix B_1 (n,inversematrix(B, n));

        Matrix C = B_1 * AA2;
        typenam matr_norm_1 = C.mnorm1();
        typenam matr_norm_Inf = C.mnormInf();
        cout << "Matrix C norm1: " << matr_norm_1 << endl;
        cout << "Matrix C norm_inf: " << matr_norm_Inf << endl;
        vector<typenam>x_k(n, 0);
        vector<typenam> x_k_x_0 = x_0;
        typenam Usum = 0;
        typenam Lsum = 0;
        //Достаточное условие 
        if (matr_norm_Inf>=1) {
            cout << "Dostatochnoe uslovie ne vypolneno" << endl;
            exit(1);
        }
        while (norm1(x_k_x_0) >= (1 - matr_norm_Inf) / matr_norm_Inf * epsilon) {
            //x_k = x_0;
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i; j < n; ++j) {
                    if (i != j) { Usum += data[i][j] * x_0[j]; }
                }                   
                x_k[i] = (b[i] - Usum) / data[i][i];
                for (size_t j = 0; j < i; ++j) {
                    if (i != j) { Lsum += data[i][j] * x_k[j]; }
                }
                x_k[i] -= Lsum / data[i][i];

                Usum = 0;
                Lsum = 0;
            }
            for (size_t k = 0; k < n; ++k) {
                x_k_x_0[k] = x_0[k] - x_k[k];
            }
            x_0 = x_k;
            ++iter_culc;
            if (iter_culc > 1000) {
                cout << "TOO MANY ITERATIONS" << endl;
                exit(1);
            }
        }
        cout << "A priori estimate iterations: " << int(log(epsilon) / log(matr_norm_Inf )) << endl;
        cout << "A posteriori estimate iterations: " << int(log(epsilon) / log(matr_norm_Inf /  (1-matr_norm_Inf))) << endl;
        cout << "Iterations: " << iter_culc << endl;
        return x_k;
    }


    

};

