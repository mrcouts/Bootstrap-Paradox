#ifndef ACCELERATION_H
#define ACCELERATION_H
#include <cmath>
#include <algorithm>
#include <armadillo>
#include "SomeUtilities.h"
#include "Dy.h"
#include "Serial.h"
#include "_5R_.h"
#include "FLControlLaw.h"

using namespace std;
using namespace arma;

class Acceleration{
public:
    Acceleration(int dof, Serial *R, FLControlLaw *u);
    Acceleration(int dof, _5R_ *R2, FLControlLaw *u);
    ~Acceleration();
    field<vec> Doit(double t, vec q0_, vec q1_);
    vec f_(double t, vec y_);
    field<vec> f2_(double t, vec y_);
    int dof;
    Serial *R;
    _5R_ *R2;
    int caso;
    FLControlLaw *u; };

#endif