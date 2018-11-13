#include <iostream>
#include <iomanip>
#include <string>
#include "SomeUtilities.h"
#include "Serial.h"
#include "RR.h"
#include "FLControlLaw.h"
#include "Acceleration.h"
#include "RK.h"
#include "Filter.h"

vec r_(double t){
    double w1 = 10;
    double w2 = 15;
    return {sin(w1*t),sin(w2*t)}; }

vec dr_(double t){
    double w1 = 10;
    double w2 = 15;
    return {w1*cos(w1*t),w2*cos(w2*t)}; }

vec d2r_(double t){
    double w1 = 10;
    double w2 = 15;
    return {-w1*w1*sin(w1*t), -w2*w2*sin(w2*t)}; }

Dy* (dy_comp)(vec q0_, vec q1_){
    Dy *dy;
    dy = new Dy(2);
    return dy; }

int main(void){
    int dof = 2;
    cube I__; I__.zeros(3,3,dof);
    I__.slice(0) << 0 << 0      << 0      << endr
                 << 0 << 0.0001 << 0      << endr
                 << 0 << 0      << 0.0001 << endr;

    I__.slice(1) << 0 << 0      << 0      << endr
                 << 0 << 0.0001 << 0      << endr
                 << 0 << 0      << 0.0001 << endr;

    Serial RR = Serial(2, {0.1, 0.1}, {0.05, 0.05},{0.1, 0.1}, I__ , {0, -9.8,0}, &fDH_RR);
    ControlLaw FL = ControlLaw(dof, 40.0, &r_, &dr_, &d2r_, &RR);
    //ControlLaw FL = ControlLaw(dof, 40.0, 40.0, 0.01, &r_, &dr_, &d2r_, &RR);

    //Acceleration AC = Acceleration(dof, &RR, &FL);
    //vec u; u.zeros(dof);
    //u = FL.Doit(0, r_(0.1), dr_(0.1));
    //cout << u << endl;
    //vec v; v.zeros(2*dof);
    //v = AC.f_(0, join_vert(r_(0.1), dr_(0.1)) );
    //cout << v << endl;

    RK rk = RK("RK8", &RR, &FL, 2);
    rk.Doit(0.001, 5.0, join_vert(r_(0.0), dr_(0.0)) );
    for(uint i = 0; i< rk.t_.n_rows; i++)
        //cout << rk.t_(i) << "; " << rk.u__(0,0,i)  << "; " << rk.u__(1,0,i) << "; " << endl;
        //cout << rk.t_(i) << "; " << r_(rk.t_(i))(1) - rk.y__(1,0,i)  << "; " << dr_(rk.t_(i))(1) - rk.y__(3,0,i) << "; " << endl;
        cout << rk.t_(i) << "; " << r_(rk.t_(i))(0) - rk.y__(0,0,i)  << "; " << r_(rk.t_(i))(1) - rk.y__(1,0,i) << "; " << endl;
        //cout << rk.t_(i) << "; " << rk.y__(0,0,i)  << "; " << rk.y__(1,0,i) << "; " << endl;

    return 0; }