#ifndef RK_H
#define RK_H
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

class RK {
public:
    RK(string method);
    ~RK();
    void Doit(double h, double tf, vec y0_, vec (*f_)(double, vec));

    int N;
    vec *a__;
    vec b_;
    vec c_;

    vec t_;
    cube y__;
    cube u__;
private:
    cube k__; };

#endif