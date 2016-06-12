#ifndef ACCELERATION_H
#define ACCELERATION_H
#include <cmath>
#include <algorithm>
#include <armadillo>
#include "SomeUtilities.h"
#include "Dy.h"
#include "Serial.h"
#include "FLControlLaw.h"

using namespace std;
using namespace arma;

class Acceleration{
public:
    Acceleration(int dof, Serial *R, FLControlLaw *u);
    ~Acceleration();
    vec Doit(double t, vec q0_, vec q1_);
    int dof;
    Serial *R;
    FLControlLaw *u; };

#endif