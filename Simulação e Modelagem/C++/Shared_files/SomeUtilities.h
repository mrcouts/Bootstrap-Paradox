#ifndef SOMEUTILITIES_H
#define SOMEUTILITIES_H
#include <cmath>
#include <armadillo>

#ifndef PI
#define PI 3.14159265358979323846264338327950288
#endif

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
mat join_diag(mat A, mat B);
vec join_vert(vec **lista, uint n);
mat join_vert(mat **lista, uint n);
mat join_diag(mat **lista, uint n);

#endif