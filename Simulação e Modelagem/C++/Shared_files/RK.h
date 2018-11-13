#ifndef RK_H
#define RK_H
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <string>
#include "Acceleration.h"
#include "GNR.h"
#include "Filter.h"


using namespace std;
using namespace arma;

class RK {
public:
    RK(string method, vec (*f_)(double, vec));
    RK(string method, Acceleration *AC);
    RK(string method, GNR *gnr);
    RK(string method, Serial *R, ControlLaw *CL, int nh);
    ~RK();
    void Doit(double h, double tf, vec y0_);

    int N;
    vec *a__;
    vec b_;
    vec c_;
    vec (*f_)(double, vec);
    Acceleration *AC;
    GNR *gnr;
    Serial *R;
    ControlLaw *CL;
    int nh;
    int counter;
    int caso;

    vec t_;
    cube y__;
    cube dy__;
    cube u__;
    cube z__;
    cube dz__;
    Filter *Fltr;
private:
    void SetMethod(string method);
    cube k__;
};

#endif