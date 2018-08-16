#ifndef  FLCONTROLLAW_H  
#define  FLCONTROLLAW_H
#include <cmath>
#include <algorithm>
#include <armadillo>
#include "Dy.h"
#include "Serial.h"
#include "Parallel.h"
#include "Reference.h"


using namespace std;
using namespace arma;

class ControlLaw {
public:
    ControlLaw(uint dof, double lambda, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Serial *R);
    ControlLaw(uint dof, double lambda, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2);
    ControlLaw(uint dof, double lambda, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Dy* (*dy_comp)(vec, vec));
    ControlLaw(uint dof, double lambda, Reference *RefObj, Serial *R);
    ControlLaw(uint dof, double lambda, Reference *RefObj, Parallel *R2);
    ControlLaw(uint dof, double lambda, Reference *RefObj, Dy* (*dy_comp)(vec, vec));

    ControlLaw(uint dof, double lambda, double K, double phi, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Serial *R);
    ControlLaw(uint dof, double lambda, double K, double phi, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2);
    ControlLaw(uint dof, double lambda, double K, double phi, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Dy* (*dy_comp)(vec, vec));
    ControlLaw(uint dof, double lambda, double K, double phi, Reference *RefObj, Serial *R);
    ControlLaw(uint dof, double lambda, double K, double phi, Reference *RefObj, Parallel *R2);
    ControlLaw(uint dof, double lambda, double K, double phi, Reference *RefObj, Dy* (*dy_comp)(vec, vec));

    ~ControlLaw();
    vec Doit(double t, vec q0_, vec q1_);

    uint dof;
    double lambda;
    double K;
    double phi;
    vec (*r_)(double);
    vec (*dr_)(double);
    vec (*d2r_)(double);
    Serial *R;
    Parallel *R2;
    Dy* (*dy_comp)(vec, vec);
    int caso;
    bool smc;
    Reference *RefObj;
    bool RefObjFlag;
};

class FLControlLaw {
public:
    FLControlLaw(uint dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Serial *R);
    FLControlLaw(uint dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2);
    FLControlLaw(uint dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Dy* (*dy_comp)(vec, vec));
    FLControlLaw(uint dof, double Kp, double Kv, Reference *RefObj, Serial *R);
    FLControlLaw(uint dof, double Kp, double Kv, Reference *RefObj, Parallel *R2);
    FLControlLaw(uint dof, double Kp, double Kv, Reference *RefObj, Dy* (*dy_comp)(vec, vec));
    ~FLControlLaw();
    vec Doit(double t, vec q0_, vec q1_);

    uint dof;
    double Kp;
    double Kv;
    vec (*r_)(double);
    vec (*dr_)(double);
    vec (*d2r_)(double);
    Serial *R;
    Parallel *R2;
    Dy* (*dy_comp)(vec, vec);
    int caso;
    Reference *RefObj;
    bool RefObjFlag;
};

class SMCLaw {
public:
    SMCLaw(uint dof, double Kp, double eta, double n, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2);
    SMCLaw(uint dof, double Kp, double eta, double n, Reference *RefObj, Parallel *R2);
    SMCLaw(uint dof, double Kp, double eta, mat Lambda_, vec gamma_, double n, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2);
    SMCLaw(uint dof, double Kp, double eta, mat Lambda_, vec gamma_, double n, Reference *RefObj, Parallel *R2);
    SMCLaw(uint dof, double Kp, vec eta_, cube Lambda__, mat Gamma_, double n, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2);
    SMCLaw(uint dof, double Kp, vec eta_, cube Lambda__, mat Gamma_, double n, Reference *RefObj, Parallel *R2);
    ~SMCLaw();
    vec Doit(double t, vec q0_, vec q1_);

    uint dof;
    double Kp;
    double eta;
    mat Lambda_;
    vec gamma_;
    vec eta_;
    cube Lambda__;
    mat Gamma_;
    double n;
    vec s_;
    vec sigma_;
    mat k;
    vec k_;
    vec (*r_)(double);
    vec (*dr_)(double);
    vec (*d2r_)(double);
    Parallel *R2;
    Reference *RefObj;
    bool RefObjFlag;
    int caso;
};

#endif