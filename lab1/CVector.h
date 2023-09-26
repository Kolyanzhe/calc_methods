#pragma once
using namespace std;
#include <vector>
#include <iostream>
#include "Header.h"

vector <typenam> product(const vector<vector<typenam>> data, const vector<typenam>& b, int n);// Óěíîćĺíčĺ ěŕňđčöű íŕ ńňîëáĺö

//Ńëîćĺíčĺ 
vector<typenam> operator+(const std::vector<typenam>&, const std::vector<typenam>&);

//Âű÷čňŕíčĺ âĺęňîđîâ
vector<typenam> operator-(const vector<typenam>&, const vector<typenam>&);

//Óěíîćĺíčĺ âĺęňîđŕ íŕ ÷čńëî
vector<typenam> operator*(const vector<typenam>&, typenam);

//Óěíîćĺíčĺ âĺęňîđŕ íŕ ÷čńëî
vector<typenam> operator*(typenam, const vector<typenam>&);

//Äĺëĺíčĺ âĺęňîđŕ íŕ ÷čńëî
vector<typenam> operator/(const vector<typenam>&, typenam);

//Ńęŕë˙đíîĺ ďđîčçâĺäĺíčĺ
typenam operator*(const vector<typenam>&, const vector<typenam>&);


//Ďđčáŕâčňü ę âĺęňîđó âĺęňîđ
vector<typenam> operator+=(vector<typenam>&, const vector<typenam>&);

//Âű÷ĺńňü čç âĺęňîđŕ âĺęňîđ
vector<typenam> operator-=(vector<typenam>&, const vector<typenam>&);

//Äîěíîćčňü âĺęňîđ íŕ ÷čńëî
vector<typenam> operator*=(vector<typenam>&, typenam);

//Đŕçäĺëčňü âĺęňîđ íŕ ÷čńëî
vector<typenam> operator/=(vector<typenam>&, typenam);

//Íîđěčđîâęŕ âĺęňîđŕ
vector<typenam> normalize(const vector<typenam>&);

//Ęâŕäđŕň ÷čńëŕ
typenam sqr(typenam);

//Ęóá ÷čńëŕ
typenam cube(typenam);


//Âűâîä â ďîňîę
ostream& operator<<(ostream&, const vector<typenam>&);

//Âűâîä íŕ ýęđŕí
void print(const vector<typenam>&);