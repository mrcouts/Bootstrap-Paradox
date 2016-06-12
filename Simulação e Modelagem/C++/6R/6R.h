#ifndef _6R_H
#define _6R_H
#include <cmath>
#include <algorithm>
#include <armadillo>

#define PI 3.14159265358979323846264338327950288

using namespace std;
using namespace arma;

mat fDH_6R(vec q0_, vec l_, vec lg_);

#endif