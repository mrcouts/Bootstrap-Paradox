#ifndef _5R__H
#define _5R__H

#include "Serial.h"
#include "SomeUtilities.h"

class _5R_:public Mecanismo{
public:
	_5R_(uint dof, Mecanismo *P, Serial **R_, uint nR_, Col<uint> u_nzi_, mat D_, mat E_, mat F_, vec f_);
    _5R_(uint dof, Mecanismo *P, Serial **R_, uint nR_, Col<uint> u_nzi_, mat D_, mat E_, mat F_, vec f_, mat P_, mat Q_, mat S_);
    ~_5R_();
    void Construct(uint dof, Mecanismo *P, Serial **R_, uint nR_, Col<uint> u_nzi_, mat D_, mat E_, mat F_, vec f_);
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
    mat **o__;
    mat **Jv_n__;
    mat **Jw_n__;
    vec **a_co_n__;
    vec **dw_co_n__;
    vec _q_;
    mat Ah_;
    mat Ao_;
    mat A_;
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