#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

#define typenam double
const typenam epsilon = 1e-4;

#include "Matrix.h"
#include"CVector.h"
typenam norm1(const vector<typenam>& a); //Îêòàýäðè÷åñêàÿ íîðìà
typenam norm2(const vector<typenam>& a); // Åâêëèäîâà íîðìà 
typenam normInf(const vector<typenam>& a); // Êóáè÷åñêàÿ íîðìà
vector <typenam> productmv(const vector<vector<typenam>> data, const vector<typenam>& b, int n);  // Óìíîæåíèå ìàòðèöû íà ñòîëáåö
vector<vector<typenam>> product(const vector<vector<typenam>>& a, const vector<vector<typenam>>& b, int n); //Óìíîæåíèå ìàòðèö
vector<vector<typenam>> inversematrix(const vector<vector<typenam>>& data, int n);

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
        cout << endl;
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
    void transpose() {
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                swap(data[i][j], data[j][i]);
    }

    //Ìàòðèöà + Ìàòðèöà
    Matrix operator+=(const Matrix& matr) {

        for (size_t row = 0; row < n; row++)
            for (size_t col = 0; col < n; col++)
                data[row][col] += matr.data[row][col];
        return *this;
    }

    //Ìàòðèöà - Ìàòðèöà
    Matrix operator-=(const Matrix& matr) {

        for (size_t row = 0; row < n; row++)
            for (size_t col = 0; col < n; col++)
                data[row][col] -= matr.data[row][col];
        return *this;
    }


    // Ñëîæåíèå ìàòðèö
    Matrix operator+(const Matrix& matr) const {
        Matrix result(*this);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                result.data[i][j] += matr.data[i][j];
        return result;
    }

    //Âû÷èòàíèå ìàòðèö
    Matrix operator-(const Matrix& matr) const {
        Matrix result(*this);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                result.data[i][j] -= matr.data[i][j];
        return result;
    }
    Matrix operator-() const {
        Matrix result(n, vector<vector<typenam>>(n, vector<typenam>(n, 0))); // Создаем новую матрицу такого же размера
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                result.data[i][j] = -data[i][j];
            }
        }
        return result;
    }
    //Óìíîæåíèå ìàòðèö
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
    std::vector<typenam> operator*(const std::vector<typenam>& vec) const {
        if (vec.size() != n) {
            throw std::runtime_error("Размерность вектора не соответствует размерности матрицы");
        }

        std::vector<typenam> result(n, 0);

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                result[i] += data[i][j] * vec[j];
            }
        }

        return result;
    }
    //Óìíîæåíèå ìàòðèöû íà ÷èñëî
    Matrix operator*(typenam k) const {
        Matrix result(*this);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                result.data[i][j] *= k;
        return result;
    }

    // Îêòàýäðè÷åñêàÿ ìàòðè÷íàÿ íîðìà 
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

    //Êóáè÷åñêàÿ ìàòðè÷íàÿ íîðìà
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
        // Èíèöèàëèçèðóåì ìàòðèöó Q êàê åäèíè÷íóþ ìàòðèöó
        vector< vector<typenam>>Q(n, vector<typenam>(n, 0.0));
        for (size_t i = 0; i < n; ++i) {
            Q[i][i] = 1.0;
        }
        if (!isUpperTriangular) {


            typenam a, c, s, buf;
            int ii, jj;
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
        Matrix QQ(n, Q);
        QQ.transpose();
        return pair<vector< vector<typenam>>, vector<vector<typenam>>>(R, QQ.data);
    }
    vector<typenam> qrDecomposition(vector<vector<typenam>>& QQ, vector<vector<typenam>>& RR, vector<typenam>& b) {
        //Считаем R,Q
        // pair<vector<vector<typenam>>, vector<vector<typenam>>> result = qrDec();
        vector<vector<typenam>> R = RR;
        vector<vector<typenam>> Q = QQ;
        //cout << "Matrix Q" << endl;
       //print(Q);
        //cout << "Matrix R" << endl;
        //print(R);
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
        // cout << isUpperTriangular << endl;
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
   
    vector<typenam> QR_lambda(int max_iterations = 100, double eps = epsilon, bool print10 = false) {
        Matrix A(n, data);
        double max;
        int iter;
        pair<vector< vector<typenam>>, vector<vector<typenam>>>RQ;
        for (int i = 0; i < max_iterations; ++i) {
            if (print10 == true) {
                if (i < 10) {
                    A.print();
                }
            }
            RQ = A.qrDec();
            Matrix R(n, RQ.first);
            Matrix Q(n, RQ.second);
            A = R * Q;
            max = fabs(A.data[1][0]);
            for (int i = 2; i < n; i++) {
                if (fabs(A.data[i][i - 1]) > max) {
                    max = fabs(A.data[i][i - 1]);
                }
            }
            if (max < eps) {
                cout << max << endl;
                iter = i;
                break;
            }
            
        }
        cout << "Iterations " << iter << endl;
        vector<typenam> lambdas(n, 0);
        for (size_t i = 0; i < n; ++i) {
            lambdas[i] = A.data[i][i];
        }
        return lambdas;
    }

    vector<typenam> QR_lambda_sdvig(int max_iterations = 100, double eps = epsilon, double delta = 1e-2, bool print10 = false) {
        Matrix A(n, data);
        vector<vector<typenam>> zeroMatrix(n, vector<typenam>(n, 0.0));
        for (int j = 0; j < n; ++j) {
            zeroMatrix[j][j] = A.data[n-1][n-1];
        }
        Matrix E(n, zeroMatrix);
        double max;
        int k = n;
        pair<vector< vector<typenam>>, vector<vector<typenam>>>RQ;
        int iterations = 0;
        while (k>1){

            for (int j = 0; j < n; ++j) {
                zeroMatrix[j][j] = A.data[k - 1][k - 1];
            }
            for (int i = 0; i < max_iterations; ++i) {

                if (print10 == true) {
                    if (i < 10) {
                        A.print();
                    }
                }
                A = A - E;
                RQ = A.qrDec();
                Matrix R(n, RQ.first);
                Matrix Q(n, RQ.second);
                A = R * Q + E;
                max = fabs(A.data[k - 1][0]);
                for (int l = 1; l < k - 1; l++) {
                    if (fabs(A.data[k - 1][l]) > max) {
                        max = fabs(A.data[k - 1][l]);
                    }
                }
                if (max < eps) {
                    k--; 
                    iterations += i;
                    break;
                }
                cout << i << endl;
            }  
            cout <<"Iterations " << iterations << endl;
        }
        vector<typenam> lambdas;
        for (size_t i = 0; i < n; ++i) {
            lambdas.push_back(A.data[i][i]);
        }
        return lambdas;
    }

    Matrix Hessenberging(double eps = epsilon) {
        vector < vector < double>> Eid(n, vector<double>(n, 0));
        for (size_t i = 0; i < n; ++i) { Eid[i][i] = 1; }
        Matrix E(n, Eid);
        Matrix A(n, data);
        double alpha = 0;
        double beta = 0;
        int iter = 0;
        Matrix T(n, Eid);

        for (size_t k = 1; k < n - 1; ++k) {
            for (size_t l = k + 1; l < n; ++l) {
                alpha = A.data[k][k - 1] / sqrt(A.data[k][k - 1] * A.data[k][k - 1] + A.data[l][k - 1] * A.data[l][k - 1]);
                beta = A.data[l][k - 1] / sqrt(A.data[k][k - 1] * A.data[k][k - 1] + A.data[l][k - 1] * A.data[l][k - 1]);
                // cout « "a=" « alpha « " beta=" « beta « endl;
                T.data[l][l] = alpha;
                T.data[k][k] = alpha;
                T.data[l][k] = -beta;
                T.data[k][l] = beta;
                //T.print();

                A = T * A;
                T.transpose();
                A = A * T;
                //A.print();
                T = E;
            }
        }
        return A;
    }

};
