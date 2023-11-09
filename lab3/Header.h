#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <ostream>
using namespace std;
#define typenam double
class Coord_vector {
private:
	int n;					//Кол-во узлов сетки
	typenam x_start;		//Начало интервала
	typenam x_end;			//Конец интервала
	typenam(*func)(typenam);//Функция
public:
	Coord_vector(int nn, typenam start, typenam end, typenam(*f)(typenam)) :n(nn), x_start(start), x_end(end), func(f) {};
	//Равномерная сетка
	vector<pair<typenam, typenam>> UniformGrid() {
		vector<pair<typenam, typenam>> points;
		typenam step = (x_end - x_start) / n;
		for (size_t i = 0; i <= n; ++i) {
			typenam x = x_start + i * step;
			typenam y = func(x);
			points.push_back(std::make_pair(x, y));
		}
		cout << "UniformGrid" << endl;
		for (size_t i = 0; i < points.size(); ++i) {
			cout << points[i].first << " " << points[i].second << endl;
		}
		cout << endl;
		return points;
	}
	//Чебышевская сетка
	vector<pair<typenam, typenam>> ChebyshevGrid() {
		vector<pair<typenam, typenam>> points;
		typenam step = (x_end - x_start) / n;
		for (size_t i = 0; i <= n; ++i) {
			typenam x = (x_start + x_end) / 2 + (x_end - x_start) / 2 * cos((2 * i + 1) * 3.14 / (2 * (n + 1)));
			typenam y = func(x);
			points.push_back(std::make_pair(x, y));
		}
		cout << "ChebyshevGrid" << endl;
		for (size_t i = 0; i < points.size(); ++i) {
			cout << points[i].first << " " << points[i].second << endl;
		}
		cout << endl;
		return points;
	}
	//Полином Лагранжа, на вход вектор точек, в которых хотим вычислять значения
	vector<pair<typenam, typenam>> LagrangePolinome(vector<typenam> x, int flag = 0) {
		vector<pair<typenam, typenam>> points;
		vector<pair<typenam, typenam>> polinome;
		if (flag == 0) {
			points = UniformGrid();
		}
		else {
			points = ChebyshevGrid();
		}
		for (size_t k = 0; k < x.size(); ++k) {
			typenam result = 0;
			typenam P = 0;
			for (int i = 0; i < n + 1; i++) {
				P = 1.0;
				for (int j = 0; j < n + 1; j++) {
					if (j != i) {
						P *= (x[k] - points[j].first) / (points[i].first - points[j].first);
					}

				}
				result += P * points[i].second;
			}
			polinome.push_back(make_pair(x[k], result));
		}
		for (size_t i = 0; i < x.size(); ++i) {
			cout << polinome[i].first << " ";
		}
		cout << endl;
		for (size_t i = 0; i < x.size(); ++i) {
			cout << polinome[i].second << " ";
		}
		cout << endl;
		//Считаем погрешность
		typenam max(0);
		for (size_t i = 0; i < x.size(); ++i) {
			typenam error = fabs(polinome[i].second - func(polinome[i].first));
			if (max < error) {
				max = error;
			}
		}
		cout << "Error: " << max << endl;
		return polinome;
	}



};
typenam func0(typenam x) {
	return x;
};
typenam func1(typenam x) {
	return x*x;
};

typenam func2(typenam x) {
	return 1/(1+pow(x, 2));
};

typenam func3(typenam x) {
	return 1 / atan(1+ 10 * pow(x, 2));
};

typenam func4(typenam x) {
	return pow((4 * pow(x, 3) + 2 * pow(x, 2) - 4 * x + 2), sqrt(2)) + asin(1 / (5 + x - x * x)) - 5;
};