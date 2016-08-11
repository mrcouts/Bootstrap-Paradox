#ifndef RK_H
#define RK_H
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <string>
#include "Acceleration.h"
#include "GNR.h"


using namespace std;
using namespace arma;

class RK {
public:
    RK(string method, vec (*f_)(double, vec));
    RK(string method, Acceleration *AC);
    RK(string method, GNR *gnr);
    ~RK();
    void Doit(double h, double tf, vec y0_);

    int N;
    vec *a__;
    vec b_;
    vec c_;
    vec (*f_)(double, vec);
    Acceleration *AC;
    GNR *gnr;
    int caso;

    vec t_;
    cube y__;
    cube u__;
private:
    void SetMethod(string method);
    cube k__;
};

#endif