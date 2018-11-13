#ifndef  FILTER_H  
#define  FILTER_H
#include <cmath>
#include <algorithm>
#include <armadillo>

using namespace std;
using namespace arma;

class Filter {
public:
    Filter(int size, mat AB_);
    ~Filter();
    vec Doit(vec u_);

    int order;
    int size;
    mat u__;
    mat y__;
    vec a_;
    vec b_;
};

mat Tustin(double T, double w0, mat ABs_);
mat Bessel(double w, int order);
mat dBessel(double w, int order);
mat Bessel_d(double T, double w, int order);
mat dBessel_d(double T, double w, int order);

#endif