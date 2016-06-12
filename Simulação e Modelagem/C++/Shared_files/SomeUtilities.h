#ifndef SOMEUTILITIES_H
#define SOMEUTILITIES_H
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

mat Rotx(double theta);
mat Roty(double theta);
mat Rotz(double theta);
mat Hx(double theta, double dx, double dy, double dz);
mat Hy(double theta, double dx, double dy, double dz);
mat Hz(double theta, double dx, double dy, double dz);
mat Ht(double dx, double dy, double dz);
mat H(char axis, double theta, double dx, double dy, double dz);
mat H_d(double a, double alpha, double d, double theta);

#endif