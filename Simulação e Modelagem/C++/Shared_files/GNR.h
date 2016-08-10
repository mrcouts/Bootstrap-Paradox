#ifndef GNR_H
#define GNR_H
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

class GNR {
public:
    GNR(string method, vec x0_, vec (*f_)(vec), mat (*J_)(vec), double tol, uint nmax);
    ~GNR();
    void Doit();

    string method;
    vec x0_;
    vec (*f_)(vec);
    mat (*J_)(vec);
    double tol;
    uint nmax;
    vec x_;
    bool convergiu;
    vec res_;
    uint n;
};

#endif