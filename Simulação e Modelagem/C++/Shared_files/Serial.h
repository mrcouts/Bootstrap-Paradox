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
};

class Serial:public Mecanismo {
public:
    Serial(int dof, vec l_, vec lg_, vec m_, cube I__, vec g_, mat (*fDH)(vec, vec, vec));
    ~Serial();
    Dy* Doit(vec q0_, vec q1_);

    vec l_;
    vec lg_;
    vec m_;
    cube I__;
    vec g_;
    cube Hr__;
    cube H__;
    cube z__;
    cube o__;
    cube og__;
    cube Jv__;
    cube Jw__;
    cube Jw2__;
    mat Jv_n_;
    mat Jw_n_;
    mat Mh_;
    vec vh_;
    vec gh_;
    //Dy *dy;
    cube w_rel__;
    cube w_arr__;
    cube dw_co__;
    cube dw_co2__;
    cube a_cen__;
    cube a_cor__;
    cube a_co__;
    vec dw_co_n_;
    vec a_co_n_;
    mat (*fDH)(vec, vec, vec); };

#endif