#include <iostream>
#include <iomanip>
#include <string>
#include "SomeUtilities.h"
#include "Serial.h"
#include "RR.h"
#include "FLControlLaw.h"
#include "Acceleration.h"
#include "RK.h"
#include "_5R_.h"

vec r_(double t){
    double r = 0.08;
    double x0 = 0.0;
    double y0 = 0.16;
    double w = 1/r;
    return {x0+r*cos(w*t), y0+r*sin(w*t)}; }

vec dr_(double t){
    double r = 0.08;
    double w = 1/r;
    return {-w*r*sin(w*t), w*r*cos(w*t)}; }

vec d2r_(double t){
    double r = 0.08;
    double w = 1/r;
    return {-w*w*r*cos(w*t), -w*w*r*sin(w*t)}; }

Dy* (dy_comp)(vec q0_, vec q1_){
    Dy *dy;
    dy = new Dy(2);
    return dy; }

int main(void){
    /*
    int dof = 2;
    cube I__; I__.zeros(3,3,dof);
    I__.slice(0) << 0 << 0        << 0        << endr
                 << 0 << 171.6e-6 << 0        << endr
                 << 0 << 0        << 171.6e-6 << endr;

    I__.slice(1) << 0 << 0        << 0        << endr
                 << 0 << 320.6e-6 << 0        << endr
                 << 0 << 0        << 320.6e-6 << endr;

    Serial RR  = Serial(2, {0.12, 0.15}, {0.5*0.12, 0.5*0.15},{0.143, 0.171}, I__ , {0, -9.8,0}, &fDH_RR);

    FLControlLaw FL = FLControlLaw(dof, 100.0, 20.0, &r_, &dr_, &d2r_, &RR);
    Acceleration AC = Acceleration(dof, &RR, &FL);
    //vec u; u.zeros(dof);
    //u = FL.Doit(0, r_(0.1), dr_(0.1));
    //cout << u << endl;
    //vec v; v.zeros(2*dof);
    //v = AC.f_(0, join_vert(r_(0.1), dr_(0.1)) );
    //cout << v << endl;

    RK rk = RK("RK8", &AC);
    rk.Doit(0.001, 0.01, join_vert(r_(0.0), dr_(0.0)) );
    for(uint i = 0; i< rk.t_.n_rows; i++)
        cout << rk.t_(i) << "; " << rk.u__(0,0,i)  << "; " << rk.u__(1,0,i) << "; " << endl;
        //cout << rk.t_(i) << "; " << r_(rk.t_(i))(1) - rk.y__(1,0,i)  << "; " << dr_(rk.t_(i))(1) - rk.y__(3,0,i) << "; " << endl;

    //cout << RR.o__.slice(0) << endl;

    */

    _5R_ Robot = _5R_();
    Robot.Doit({0.0, 0.20, 1.01377246756945, 1.41460649673445, 1.01377246756945, 1.41460649673445},{0,0,0,0,0,0});

    //cout << Robot.dy->Mh_ << endl;
    //cout << Robot.dy->vh_ << endl;
    //cout << Robot.dy->gh_ << endl;

    FLControlLaw FL = FLControlLaw(6, 100.0, 20.0, &r_, &dr_, &d2r_, &Robot);
    Acceleration AC = Acceleration(6, &Robot, &FL);

    //vec u; u.zeros(2);
    //u = FL.Doit(0, {0.08,  0.16, 0.305030291698133, 1.86386236511897, 1.45111035931733, 1.41460649673445}, zeros(6));
    //cout << u << endl;

    vec x0_ = {0.08,  0.16, 0.305030291698133, 1.86386236511897, 1.45111035931733, 1.41460649673445, 0, 0, 0, 0, 0, 0};

    //vec v; v.zeros(2*6);
    //v = AC.f_(0, x0_ );
    //cout << v << endl;

    //field<vec> F(2);
    //F = AC.Doit(0, {0.08,  0.16, 0.305030291698133, 1.86386236511897, 1.45111035931733, 1.41460649673445}, {0, 0, 0, 0, 0, 0} );
    //cout << F << endl;

    RK rk = RK("RK8", &AC);
    rk.Doit(0.001, 4*1.2, x0_);
    for(uint i = 0; i< rk.t_.n_rows; i++)
        cout << rk.t_(i) << "; " << rk.u__(0,0,i)  << "; " << rk.u__(1,0,i) << "; " << endl;

    return 0;
}