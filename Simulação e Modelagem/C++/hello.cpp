#include <iostream>
#include <iomanip>
#include <string>
#include "SomeUtilities.h"
#include "Serial.h"
#include "RR.h"
#include "6R.h"
#include "RK.h"

//Simple C++ program
int main(void)
{
	//cout << setprecision(16);
    cout << "Hello! This is a C++ program." << endl;
    cout << H_d(1, 0.1, 2, 0.2) << endl;

    /*RR
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
    	t += 0.01;
    }
    */

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

    //bool array[] = {true, true};
    Serial R6 = Serial(dof, l_, lg_ , m_, I__ , {0, -9.8,0}, &fDH_6R);

    double t = 0;
    vec q0_; q0_.zeros(dof);
    vec q1_; q1_.zeros(dof);
    vec q2_; q2_.zeros(dof);
    double w1 = 10;
    double w2 = 15;
    for(int i = 0; i<2; i++){
    	q0_ = {sin(w1*t),sin(w2*t),sin(w1*t),sin(w2*t),sin(w1*t),sin(w2*t)};
    	q1_ = {w1*cos(w1*t),w2*cos(w2*t),w1*cos(w1*t),w2*cos(w2*t),w1*cos(w1*t),w2*cos(w2*t)};
    	q2_ = {-w1*w1*sin(w1*t), -w2*w2*sin(w2*t),-w1*w1*sin(w1*t), -w2*w2*sin(w2*t),-w1*w1*sin(w1*t), -w2*w2*sin(w2*t)};
    	R6.Doit(q0_,q1_);
    	cout << R6.Mh_*q2_ + R6.vh_ + R6.gh_ << endl;
    	t += 0.01;
    }

    cout << R6.dw_co_n_ << endl;
    cout << R6.a_co_n_ << endl;

    I__.clear();

    vec *vetor;
    vetor = new vec[2];
    vetor[0] = {1};
    vetor[1] = {1,2};

    cout << vetor[0] << endl;
    cout << vetor[1] << endl;

    string str1 = "RK4";
    cout << str1 << endl;

    RK RK4 = RK("RK4");
    cout << RK4.c_ << endl;

    //bool array[] = {true, true};

    return 0;
}