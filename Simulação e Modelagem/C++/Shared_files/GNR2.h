#ifndef GNR2_H
#define GNR2_H
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <string>
#include "RK.h"
#include "GNR.h"

using namespace std;
using namespace arma;

class GNR2:public GNR{
public:
    GNR2(string method, vec x0_, vec (*f_)(vec), mat (*J_)(vec), double tol, uint nmax);
    ~GNR2();
    void Doit();
    RK *rk;
    //GNR *gnr;
};

#endif