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

vector<vector<double>> R(vector<vector<typenam>> A, int n, vector<typenam> b);
vector<typenam> qrDecomposition(std::vector<std::vector<typenam>>& A, std::vector<std::vector<typenam>>& R, std::vector<std::vector<typenam>>& Q, vector<typenam> b);