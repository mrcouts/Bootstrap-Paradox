#ifndef  SMCLAW_H  
#define  SMCLAW_H
#include <cmath>
#include <algorithm>
#include <armadillo>
#include "Dy.h"
#include "Serial.h"
#include "Parallel.h"
#include "Reference.h"


using namespace std;
using namespace arma;

class SMCLaw {
public:
    SMCLaw(int dof, double Kp, double eta, mat K_, vec k_, double n, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2);
    SMCLaw(int dof, double Kp, double eta, mat K_, vec k_, double n, Reference *RefObj, Parallel *R2);
    ~SMCLaw();
    vec Doit(double t, vec q0_, vec q1_);

    int dof;
    double Kp;
    double eta;
    mat K_;
    vec k_;
    double n;
    vec s_;
    vec sigma_;
    mat k;
    vec (*r_)(double);
    vec (*dr_)(double);
    vec (*d2r_)(double);
    Parallel *R2;
    Reference *RefObj;
    bool RefObjFlag;
};

#endif