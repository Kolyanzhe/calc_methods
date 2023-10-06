#pragma once
#pragma once
using namespace std;
#include <vector>
#include <iostream>
#include "Header.h"

vector <typenam> product(const vector<vector<typenam>> data, const vector<typenam>& b, int n);// Умножение матрицы на столбец

//Сложение 
vector<typenam> operator+(const std::vector<typenam>&, const std::vector<typenam>&);

//Вычитание векторов
vector<typenam> operator-(const vector<typenam>& a1, const vector<typenam>& a2);

//Умножение вектора на число
vector<typenam> operator*(const vector<typenam>&, typenam);

//Умножение вектора на число
vector<typenam> operator*(typenam, const vector<typenam>&);

//Деление вектора на число
vector<typenam> operator/(const vector<typenam>&, typenam);

//Скалярное произведение
typenam operator*(const vector<typenam>&, const vector<typenam>&);


//Прибавить к вектору вектор
vector<typenam> operator+=(vector<typenam>&, const vector<typenam>&);

//Вычесть из вектора вектор
vector<typenam> operator-=(vector<typenam>&, const vector<typenam>&);

//Домножить вектор на число
vector<typenam> operator*=(vector<typenam>&, typenam);

//Разделить вектор на число
vector<typenam> operator/=(vector<typenam>&, typenam);

//Нормировка вектора
vector<typenam> normalize(const vector<typenam>&);

//Квадрат числа
typenam sqr(typenam);

//Куб числа
typenam cube(typenam);


//Вывод в поток
ostream& operator<<(ostream&, const vector<typenam>&);

//Вывод на экран
void print(const vector<typenam>&);
