#ifndef Parallel_H
#define Parallel_H

#include "Serial.h"
#include "SomeUtilities.h"

class Parallel:public Mecanismo{
public:
	Parallel(uint dof, Mecanismo *P, Serial **R_, uint nR_, Col<uint> u_nzi_, mat D_, mat E_, mat F_, vec f_);
    Parallel(uint dof, Mecanismo *P, Serial **R_, uint nR_, Col<uint> u_nzi_, mat D_, mat E_, mat F_, vec f_, mat P_, mat Q_, mat S_);
    ~Parallel();
    void Construct(uint dof, Mecanismo *P, Serial **R_, uint nR_, Col<uint> u_nzi_, mat D_, mat E_, mat F_, vec f_);
    void DoitKin(vec q0_, vec q1_);
    Dy* Doit(vec q0_, vec q1_);

	int nq;
	mat M_;
    vec v_;
    vec g_;
    vec **q0__;
    vec **q1__;
    mat **M__;
    vec **v__;
    vec **g__;
    vec **o__;
    mat **Jv_n__;
    mat **Jw_n__;
    vec **a_co_n__;
    vec **dw_co_n__;
    vec _q_;
    mat Ah_;
    mat Ao1_;
    mat Ao2_;
    mat Ao_;
    mat A_;
    vec b1_;
    vec b2_;
    vec b_;
    mat C_;
    mat Z_;
    Mecanismo *P;
    Serial **R_;
    Col<uint> u_nzi_;
    uint nR_;
    int caso;

    mat D_;
    mat E_;
    mat F_;
    vec f_;
    mat P_;
    mat Q_;
    mat S_;
};

#endif