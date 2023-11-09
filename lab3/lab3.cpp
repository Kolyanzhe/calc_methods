#include "Header.h"

using namespace std;

int main() {
    int n = 2;
    //int n = 10;
    //int n = 30;
    Coord_vector coords(n, -1, 1, func4);
    //Просто пример
    vector<typenam> x = { -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.3, 0.1, 0.3, 0.5, 0.7, 0.9 };
    coords.LagrangePolinome(x);



    return 0;
}