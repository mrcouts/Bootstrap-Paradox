#ifndef  REFERENCE_H  
#define  REFERENCE_H
#include <cmath>
#include <algorithm>
#include <armadillo>

using namespace std;
using namespace arma;

class Reference{
public:
    Reference(double tf, vec x0_, vec xf_);
    ~Reference();
    void Doit(double t);

    uint dof;
    vec r_;
    vec dr_;
    vec d2r_;
    double t;
    double tf;
    vec x0_;
    vec dx_;
};

#endif