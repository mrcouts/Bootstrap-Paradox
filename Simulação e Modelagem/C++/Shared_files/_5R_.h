#ifndef _5R__H
#define _5R__H

#include "Serial.h"
#include "SomeUtilities.h"
#include "RR.h"

class _5R_:public Mecanismo{
public:
	_5R_();
    ~_5R_();
    Dy* Doit(vec q0_, vec q1_);

	int dof;
	int nq;
	mat M_;
    vec v_;
    vec g_;
    vec _q_;
    mat Ah_;
    mat Ao_;
    mat A_;
    vec b_;
    mat C_;
    mat Z_;
    mat Mh_;
    vec vh_;
    vec gh_;
    mat H_;
    mat h_;
    Mecanismo *P;
    Serial *RR1;
    Serial *RR2;
};

#endif