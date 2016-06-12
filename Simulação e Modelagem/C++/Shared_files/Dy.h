#ifndef DY_H
#define DY_H
#include <cmath>
#include <algorithm>
#include <armadillo>
#include "SomeUtilities.h"

using namespace std;
using namespace arma;

class Dy {
public:
    Dy(int dof);
    ~Dy();

    mat Mh_;
    vec vh_;
    vec gh_; };

#endif