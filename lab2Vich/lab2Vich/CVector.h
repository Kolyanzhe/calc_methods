#pragma once
#pragma once
using namespace std;
#include <vector>
#include <iostream>
#include "Header.h"

vector <typenam> product(const vector<vector<typenam>> data, const vector<typenam>& b, int n);// ��������� ������� �� �������

//�������� 
vector<typenam> operator+(const std::vector<typenam>&, const std::vector<typenam>&);

//��������� ��������
vector<typenam> operator-(const vector<typenam>& a1, const vector<typenam>& a2);

//��������� ������� �� �����
vector<typenam> operator*(const vector<typenam>&, typenam);

//��������� ������� �� �����
vector<typenam> operator*(typenam, const vector<typenam>&);

//������� ������� �� �����
vector<typenam> operator/(const vector<typenam>&, typenam);

//��������� ������������
typenam operator*(const vector<typenam>&, const vector<typenam>&);


//��������� � ������� ������
vector<typenam> operator+=(vector<typenam>&, const vector<typenam>&);

//������� �� ������� ������
vector<typenam> operator-=(vector<typenam>&, const vector<typenam>&);

//��������� ������ �� �����
vector<typenam> operator*=(vector<typenam>&, typenam);

//��������� ������ �� �����
vector<typenam> operator/=(vector<typenam>&, typenam);

//���������� �������
vector<typenam> normalize(const vector<typenam>&);

//������� �����
typenam sqr(typenam);

//��� �����
typenam cube(typenam);


//����� � �����
ostream& operator<<(ostream&, const vector<typenam>&);

//����� �� �����
void print(const vector<typenam>&);
