#ifndef SERIAL_H
#define SERIAL_H
#include <cmath>
#include <algorithm>
#include <armadillo>
#include "SomeUtilities.h"
#include "Dy.h"

using namespace std;
using namespace arma;

class Mecanismo{
public:
    Mecanismo(uint dof);
    ~Mecanismo();
    Dy* Doit(vec q0_, vec q1_);
    uint dof;
    Dy *dy;
    vec q0_;
    vec q1_; };

class Serial:public Mecanismo {
public:
    Serial(int dof, vec l_, vec lg_, vec m_, cube I__, vec g_, mat (*fDH)(vec, vec, vec));
    ~Serial();
    Dy* Doit(vec q0_, vec q1_);

    vec l_;
    vec lg_;
    vec m_;
    cube I__;
    cube Ig__;
    vec g_;
    cube Hr__;
    cube H__;
    mat z__;
    mat o__;
    vec x_;
    mat og__;
    cube Jv__;
    cube Jw__;
    mat Jv_n_;
    mat Jw_n_;
    mat Mh_;
    vec vh_;
    vec gh_;
    //Dy *dy;
    mat v__;
    cube v___;
    mat w__;
    vec a_co_n_;
    vec dw_co_n_;
    mat a_co_n_i_;
    mat dw_co_n_i_;
    cube a_co_ij__;
    mat a_co__;
    mat dw_co__;
    mat (*fDH)(vec, vec, vec); };

#endif