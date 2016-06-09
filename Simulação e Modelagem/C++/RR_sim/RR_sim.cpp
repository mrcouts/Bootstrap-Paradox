#include <iostream>
#include <iomanip>
#include <string>
#include "SomeUtilities.h"
#include "Serial.h"
#include "RR.h"

int main(void){
    int dof = 2;
    cube I__; I__.zeros(3,3,2);
    I__.slice(0) << 0 << 0      << 0      << endr
                 << 0 << 0.0001 << 0      << endr
                 << 0 << 0      << 0.0001 << endr;

    I__.slice(1) << 0 << 0      << 0      << endr
                 << 0 << 0.0001 << 0      << endr
                 << 0 << 0      << 0.0001 << endr;

    //bool array[] = {true, true};
    Serial RR = Serial(2, {0.1, 0.1}, {0.05, 0.05},{0.1, 0.1}, I__ , {0, -9.8,0}, &fDH_RR);

    double t = 0;
    vec q0_;
    vec q1_;
    vec q2_;
    double w1 = 10;
    double w2 = 15;
    vec u; u.zeros(dof);
    for(int i = 0; i<1000; i++){
    	q0_ = {sin(w1*t),sin(w2*t)};
    	q1_ = {w1*cos(w1*t),w2*cos(w2*t)};
    	q2_ = {-w1*w1*sin(w1*t), -w2*w2*sin(w2*t)};
    	RR.Doit(q0_,q1_);
    	u = RR.Mh_*q2_ + RR.vh_ + RR.gh_;
        cout << t << "; ";
        for(int j = 0; j<dof; j++)
            cout << u(j) << "; ";
        cout << endl;
    	t += 0.01; }

    return 0; }