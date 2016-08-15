#ifndef GNR_H
#define GNR_H
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <string>
#include "Parallel.h"

using namespace std;
using namespace arma;

class GNR {
public:
    GNR(vec (*f_)(vec), mat (*J_)(vec), double tol, uint nmax);
    GNR(Parallel *R, double tol, uint nmax);
    ~GNR();
    vec g_(vec y_);
    void Doit(vec x0_);

    Parallel *R;
    vec x0_;
    vec (*f_)(vec);
    mat (*J_)(vec);
    double tol;
    uint nmax;
    vec x_;
    bool convergiu;
    vec res_;
    uint n;
    uint caso;
};

#endif