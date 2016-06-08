#include <iostream>
#include <iomanip>
#include <string>
#include "SomeUtilities.h"
#include "Serial.h"
#include "RR.h"

int main(void){
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
    for(int i = 0; i<10000; i++){
    	q0_ = {sin(w1*t),sin(w2*t)};
    	q1_ = {w1*cos(w1*t),w2*cos(w2*t)};
    	q2_ = {-w1*w1*sin(w1*t), -w2*w2*sin(w2*t)};
    	RR.Doit(q0_,q1_);
    	cout << RR.Mh_*q2_ + RR.vh_ + RR.gh_ << endl;
    	t += 0.01; }

    return 0; }