#include <iostream>
#include <iomanip>
#include <string>
#include "SomeUtilities.h"
#include "Serial.h"
#include "PRRP.h"

int main(void) {
    int dof = 4;
    vec l_ = {0.000, 0.380, 0.000, 0.200};
    vec lg_= {0.000, 0.191, 0.000, 0.100};
    vec m_ = {0.738, 0.326, 0.136, 0.128};

    cube I__; I__.zeros(3,3,dof);
    I__.slice(0) << 0  << 0  << 0 << endr
                 << 0  << 0  << 0 << endr
                 << 0  << 0  << 0 << endr;

    I__.slice(1) << 0 << 0 << 0        << endr
                 << 0 << 0 << 0        << endr
                 << 0 << 0 << 13662e-6 << endr;

    I__.slice(2) << 0  << 0  << 0 << endr
                 << 0  << 0  << 0 << endr
                 << 0  << 0  << 0 << endr;

    I__.slice(3) << 0  << 0  << 0 << endr
                 << 0  << 0  << 0 << endr
                 << 0  << 0  << 0 << endr;

    Serial PRRP = Serial(dof, l_, lg_ , m_, I__ , {0, -9.8,0}, &fDH_PRRP);

    double t = 0;
    vec q0_; q0_.zeros(dof);
    vec q1_; q1_.zeros(dof);
    vec q2_; q2_.zeros(dof);
    double w1 = 10;
    double w2 = 15;
    vec u; u.zeros(dof);
    for(int i = 0; i<1000; i++){
    	q0_ = {sin(w1*t),sin(w2*t),sin(w1*t),sin(w2*t)};
    	q1_ = {w1*cos(w1*t),w2*cos(w2*t),w1*cos(w1*t),w2*cos(w2*t)};
    	q2_ = {-w1*w1*sin(w1*t), -w2*w2*sin(w2*t),-w1*w1*sin(w1*t), -w2*w2*sin(w2*t)};
    	PRRP.Doit(q0_,q1_);
    	u = PRRP.Mh_*q2_ + PRRP.vh_ + PRRP.gh_;
        cout << t << "; ";
        for(int j = 0; j<dof; j++)
            cout << u(j) << "; ";
        cout << endl;
    	t += 0.01;
    }

    l_.clear();
    lg_.clear();
    m_.clear();
    I__.clear();

    return 0; }