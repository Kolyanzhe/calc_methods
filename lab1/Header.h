#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;
#define typenam double

class Matrix {
private:
    int rows;
    int cols;
    vector<vector<typenam>> data;

public:
    Matrix(int rows, int cols, vector<vector<typenam>> data) {
        this->cols = cols;
        this->rows = rows;
        this->data = data;
    }

    void print() const {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                cout << data[i][j] << "    ";
            }
            cout << endl;
        }
    }
    int numRows() const { return rows; }

    void swapRows(int row1, int row2) {
        if (row1 < 0 || row1 >= rows || row2 < 0 || row2 >= rows) {
            cout << "index_error" << endl;
            exit(1);
        }
        swap(data[row1], data[row2]);
    }

    void print(vector<vector<typenam>>& matrix) {
        for (const auto& row : matrix) {
            for (const typenam& element : row) {
                std::cout << element << " ";
            }
            std::cout << std::endl;
        }
    }
    vector<typenam> qrDecomposition(vector<typenam> b) {

        vector<vector<typenam>> R(data);
        size_t n = rows;
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
        if (!isUpperTriangular) {
            // Инициализируем матрицу Q как единичную матрицу
            vector< vector<typenam>>Q(n, vector<typenam>(n, 0.0));
            for (size_t i = 0; i < n; ++i) {
                Q[i][i] = 1.0;
            }
            typenam a, c, s, r, temp;
            for (size_t j = 0; j < n; ++j) {
                for (size_t i = n - 1; i >= j + 1; --i) { // Изменено условие
                   //Во временных переменных нет ничего плохого. 
                    a = R[i - 1][j];
                    typenam b = R[i][j];
                    r = sqrt(a * a + b * b);

                    // Вычисляем c и s
                    c = a / r;
                    s = -b / r;

                    // Обновляем матрицу R
                    for (size_t k = 0; k < n; ++k) {
                        temp = c * R[i - 1][k] - s * R[i][k];
                        R[i][k] = s * R[i - 1][k] + c * R[i][k];
                        R[i - 1][k] = temp;

                    }
                    // Обновляем матрицу Q
                    for (size_t k = 0; k < n; ++k) {
                        temp = c * Q[i - 1][k] - s * Q[i][k];
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
            print(R);
            vector <typenam> x(n);
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
            print(R);
            vector <typenam> x(n);
            for (int i = n - 1; i >= 0; --i) {
                x[i] = b[i];
                for (int j = i + 1; j < n; ++j) {
                    x[i] -= R[i][j] * x[j];
                }
                x[i] /= R[i][i];
            }
        }

    }

    vector<typenam> solveGaussPartialPivoting(const vector<typenam>& b) {
        int n = numRows();

        vector<vector<typenam>> A_copy(data);
        vector<typenam> result(b);
        for (int i = 0; i < n; ++i) {
            int maxRow = i;
            for (int j = i + 1; j < n; ++j) {
                if (fabs(A_copy[j][i]) > fabs(A_copy[maxRow][i])) {
                    maxRow = j;
                }
            }

            if (i != maxRow) {
                swap(A_copy[i], A_copy[maxRow]);
                swap(result[i], result[maxRow]);
            }
       
            if (abs(A_copy[i][i]) < 1e-9) {
                cout << "Systema nesovmestna)" << endl;
                exit(1);
            }

            for (int j = i + 1; j < n; ++j) {
                typenam factor = A_copy[j][i] / A_copy[i][i];
                for (int k = i; k < n; ++k) {
                    A_copy[j][k] -= factor * A_copy[i][k];
                }
                result[j] -= factor * result[i];
            }
        }

        vector <typenam> x(n);
        for (int i = n - 1; i >= 0; --i) {
            x[i] = result[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= A_copy[i][j] * x[j];
            }
            x[i] /= A_copy[i][i];
        }
        cout << "upper triangular matrix: " << endl;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                cout << A_copy[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
        return x;
    }
};

typenam norm1(const vector<typenam>& a); //Îêòàýäðè÷åñêàÿ íîðìà
typenam norm2(const vector<typenam>& a); // Åâêëèäîâà íîðìà 
typenam normInf(const vector<typenam>& a); // Êóáè÷åñêàÿ íîðìà 
void print(vector<vector<typenam>>& matrix);

