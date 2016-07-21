#ifndef RR_H
#define RR_H
#include <cmath>
#include <algorithm>
#include <armadillo>

#ifndef PI
#define PI 3.14159265358979323846264338327950288
#endif

using namespace std;
using namespace arma;

mat fDH_RR(vec q0_, vec l_, vec lg_);

#endif