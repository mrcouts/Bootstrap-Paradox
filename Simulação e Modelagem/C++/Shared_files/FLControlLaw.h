#ifndef  FLCONTROLLAW_H  
#define  FLCONTROLLAW_H
#include <cmath>
#include <algorithm>
#include <armadillo>
#include "Dy.h"
#include "Serial.h"
#include "Parallel.h"

using namespace std;
using namespace arma;

class FLControlLaw {
public:
    FLControlLaw(int dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Serial *R);
    FLControlLaw(int dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2);
    FLControlLaw(int dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Dy* (*dy_comp)(vec, vec));
    ~FLControlLaw();
    vec Doit(double t, vec q0_, vec q1_);

    int dof;
    double Kp;
    double Kv;
    vec (*r_)(double);
    vec (*dr_)(double);
    vec (*d2r_)(double);
    Serial *R;
    Parallel *R2;
    Dy* (*dy_comp)(vec, vec);
    int caso; };

#endif