#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

class Matrix {
private:
    int rows;
    int cols;
    vector<vector<double>> data;

public:
    Matrix(int rows, int cols, vector<vector<double>> data) {
        this->cols = cols;
        this->rows = rows;
        this->data = data;
    }

    void print() const {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                cout << data[i][j] << " ";
            }
            cout << endl;
        }
    }
    int numRows() const {return rows;}

    void swapRows(int row1, int row2) {
        if (row1 < 0 || row1 >= rows || row2 < 0 || row2 >= rows) {
            cout << "index_error" << endl;
            exit(1);
        }
        swap(data[row1], data[row2]);
    }

    vector<double> solveGaussPartialPivoting(const vector<double>& b) {
        int n = numRows();

        vector<vector<double>> A_copy(data);
        vector<double> result(b);

        for (int i = 0; i < n; ++i) {
            int maxRow = i;
            for (int j = i + 1; j < n; ++j) {
                if (abs(A_copy[j][i]) > abs(A_copy[maxRow][i])) {   
                    maxRow = j;
                }
            }
            if (i != maxRow) {
                swapRows(i, maxRow);
            }

            if (abs(A_copy[i][i]) < 1e-9) {
                cout << "Matrica nesovmestna)" << endl;
                exit(1);
            }

            for (int j = i + 1; j < n; ++j) {
                double factor = A_copy[j][i] / A_copy[i][i];
                for (int k = i; k < n; ++k) {
                    A_copy[j][k] -= factor * A_copy[i][k];
                }
                result[j] -= factor * result[i];
            }
        }

        vector<double> x(n);
        for (int i = n - 1; i >= 0; --i) {
            x[i] = result[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= A_copy[i][j] * x[j];
            }
            x[i] /= A_copy[i][i];
        }
        /*for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                cout << A_copy[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;*/
        return x;
    }
};


int main() {
    int rows = 4;
    int cols = 4;
    vector<vector<double>> data1 = {
        {1,1,1,1},
        {0,1,1,1},
        {0,0,1,1},
        {0,0,0,1}
    };
    vector<vector<double>> data2 = {
    {0,0,0,1},
    {0,0,1,1},
    {0,1,1,1},
    {1,1,1,1}
    };
    vector<vector<double>> data3 = {
    {1,1,1,1},
    {2,3,3,3},
    {2,4,4,4},
    {4,5,6,7}
    };
    vector<vector<double>> data4 = {
    {10,6,2,0},
    {5,1,-2,4},
    {3,5,1,-1},
    {0,6,-2,2}
    };
    vector<vector<double>> data5 = {
    {28.859, -0.008, 2.406, 19.240},
    {14.436, -0.001, 1.203, 9.624},
    {120.204, -0.032, 10.024, 80.144},
    {-57.714, 0.016, -4.812, -38.478}
    };
    vector<double> b1({4,3,2,1});
    vector<double> b2({1,2,3,4});
    vector<double> b3({4,11,15,22});
    vector<double> b4({25,14,10,8});
    vector<double> b5({30.459, 18.248, 128.156, -60.908});
    Matrix A(rows, cols, data5);
    A.print();
    vector<double> x = A.solveGaussPartialPivoting(b5);
    A.print();
    cout << "Solution:" << endl;
    for (int i = 0; i < 4; ++i) {
        cout << "x[" << i+1 << "] = " << x[i] << endl;
    }
    return 0;
}


















/*
vector<double> gaussPartialPivoting(const Matrix& A, const vector<double>& b) {
    int n = A.numRows();

    // Создаем копии матрицы и вектора b
    vector<vector<double>> A_copy(n, vector<double>(n));
    vector<double> result(b);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A_copy[i][j] = A.get(i, j);
        }
    }

    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int j = i + 1; j < n; ++j) {
            if (abs(A_copy[j][i]) > abs(A_copy[maxRow][i])) {
                maxRow = j;
            }
        }

        if (abs(A_copy[i][i]) < EPSILON) {
            cerr << "Система уравнений вырождена или бесконечное количество решений." << endl;
            exit(1);
        }

        for (int j = i + 1; j < n; ++j) {
            double factor = A_copy[j][i] / A_copy[i][i];
            for (int k = i; k < n; ++k) {
                A_copy[j][k] -= factor * A_copy[i][k];
            }
            result[j] -= factor * result[i];
        }
    }

    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = result[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A_copy[i][j] * x[j];
        }
        x[i] /= A_copy[i][i];
    }

    return x;
    */