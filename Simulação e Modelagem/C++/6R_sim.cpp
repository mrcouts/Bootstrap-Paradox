#include <iostream>
#include <iomanip>
#include <string>
#include "SomeUtilities.h"
#include "Serial.h"
#include "RR.h"
#include "6R.h"
#include "RK.h"

int main(void) {
    int dof = 6;
    vec l_ = {0.305, 0.045, 0.045, 0.372, 0.0125, 0.010};
    vec lg_= {0.158, 0.029, 0.016, 0.186, 0.0045, 0.005};
    vec m_ = {0.100, 0.033, 0.034, 0.060, 0.023 , 0.020};

    cube I__; I__.zeros(3,3,dof);
    I__.slice(0) << 0 << 0       << 0      << endr
                 << 0 << 1165e-6 << 0      << endr
                 << 0 << 0       << 0      << endr;

    I__.slice(1) << 8.4e-6 << 0      << 0      << endr
                 << 0      << 0.8e-6 << 0      << endr
                 << 0      << 0      << 8.4e-6 << endr;

    I__.slice(2) << 8.3e-6 << 0      << 0      << endr
                 << 0      << 8.3e-6 << 0      << endr
                 << 0      << 0      << 0.8e-6 << endr;

    I__.slice(3) << 673e-6 << 0      << 0      << endr
                 << 0      << 0.8e-6 << 0      << endr
                 << 0      << 0      << 673e-6 << endr;

    I__.slice(4) << 1.3e-6 << 0      << 0      << endr
                 << 0      << 1.3e-6 << 0      << endr
                 << 0      << 0      << 0.8e-6 << endr;

    I__.slice(5) << 1.3e-6 << 0      << 0      << endr
                 << 0      << 1.3e-6 << 0      << endr
                 << 0      << 0      << 0.8e-6 << endr;

    Serial R6 = Serial(dof, l_, lg_ , m_, I__ , {0, -9.8,0}, &fDH_6R);

    double t = 0;
    vec q0_; q0_.zeros(dof);
    vec q1_; q1_.zeros(dof);
    vec q2_; q2_.zeros(dof);
    double w1 = 10;
    double w2 = 15;
    for(int i = 0; i<10; i++){
    	q0_ = {sin(w1*t),sin(w2*t),sin(w1*t),sin(w2*t),sin(w1*t),sin(w2*t)};
    	q1_ = {w1*cos(w1*t),w2*cos(w2*t),w1*cos(w1*t),w2*cos(w2*t),w1*cos(w1*t),w2*cos(w2*t)};
    	q2_ = {-w1*w1*sin(w1*t), -w2*w2*sin(w2*t),-w1*w1*sin(w1*t), -w2*w2*sin(w2*t),-w1*w1*sin(w1*t), -w2*w2*sin(w2*t)};
    	R6.Doit(q0_,q1_);
    	cout << R6.Mh_*q2_ + R6.vh_ + R6.gh_ << endl;
    	t += 0.01;
    }

    l_.clear();
    lg_.clear();
    m_.clear();
    I__.clear();

    return 0; }