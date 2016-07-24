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

    _5R_ Robot = _5R_();

    FLControlLaw FL = FLControlLaw(6, 400.0, 40.0, &r_, &dr_, &d2r_, &Robot);
    Acceleration AC = Acceleration(6, &Robot, &FL);

    vec x0_ = {0.08,  0.16, 0.305030291698133, 1.86386236511897, 1.45111035931733, 1.41460649673445, 0, 0, 0, 0, 0, 0};

    RK rk = RK("RK8", &AC);
    rk.Doit(0.001, 4*1.2, x0_);
    for(uint i = 0; i< rk.t_.n_rows; i++)
        //cout << rk.t_(i) << "; " << r_(rk.t_(i))(0) - rk.y__(0,0,i)  << "; " << r_(rk.t_(i))(1) - rk.y__(1,0,i) << "; " << endl;
        cout << rk.t_(i) << "; " << rk.u__(0,0,i)  << "; " << rk.u__(1,0,i) << "; " << endl;

    return 0;
}