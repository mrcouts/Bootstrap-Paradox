#include <iostream>
#include <iomanip>
#include <string>
#include "SomeUtilities.h"
#include "Serial.h"
#include "6R.h"
#include "PRRP.h"
#include "FLControlLaw.h"
#include "Acceleration.h"
#include "RK.h"
#include "Parallel.h"
#include "Reference.h"

vec r_(double t){
    double r = 0.08;
    double x0 = 0.0;
    double y0 = 0.16;
    double w = 1/r;
    return {x0+r*cos(w*t), y0+r*sin(w*t),0.480};
    //return {0,0,0.480};
}

vec dr_(double t){
    double r = 0.08;
    double w = 1/r;
    return {-w*r*sin(w*t), w*r*cos(w*t),0};
    //return {0,0,0};
}

vec d2r_(double t){
    double r = 0.08;
    double w = 1/r;
    return {-w*w*r*cos(w*t), -w*w*r*sin(w*t),0};
    //return {0,0,0};
}

int main(){

    Mecanismo P = Mecanismo(3);

    vec l_ = {0.305, 0.045, 0.045 , 0.372, 0.0125, 0.000};
    vec lg_= {0.158, 0.025, 0.0116, 0.175, 0.0042, 0.000};
    vec m_ = {0.100, 0.033, 0.034 , 0.060, 0.023 , 0.000};

    cube I__; I__.zeros(3,3,6);
    I__.slice(0) << 1154e-6 << 0       << 0       << endr
                 << 0       << 1165e-6 << 0       << endr
                 << 0       << 0       << 12.7e-6 << endr;

    I__.slice(1) << 8.0e-6 << 0      << 0      << endr
                 << 0      << 1.5e-6 << 0      << endr
                 << 0      << 0      << 8.2e-6 << endr;

    I__.slice(2) << 7.9e-6 << 0      << 0      << endr
                 << 0      << 8.2e-6 << 0      << endr
                 << 0      << 0      << 1.2e-6 << endr;

    I__.slice(3) << 672.8e-6 << 0      << 0        << endr
                 << 0        << 1.7e-6 << 0        << endr
                 << 0        << 0      << 672.8e-6 << endr;

    I__.slice(4) << 1.6e-6 << 0      << 0      << endr
                 << 0      << 1.3e-6 << 0      << endr
                 << 0      << 0      << 1.2e-6 << endr;

    I__.slice(5) << 0 << 0 << 0 << endr
                 << 0 << 0 << 0 << endr
                 << 0 << 0 << 0 << endr;

    Serial _6R1 = Serial(6, l_, lg_ , m_, I__ , {0, 9.8,0}, &fDH_6R);
    Serial _6R2 = Serial(6, l_, lg_ , m_, I__ , {0, 9.8,0}, &fDH_6R);

    vec l2_ = {0.000, 0.380, 0.000, 0.200};
    vec lg2_= {0.000, 0.191, 0.000, 0.100};
    vec m2_ = {0.738, 0.326, 0.136, 0.128};

    cube I2__; I2__.zeros(3,3,4);
    I2__.slice(0) << 0  << 0  << 0 << endr
                  << 0  << 0  << 0 << endr
                  << 0  << 0  << 0 << endr;

    I2__.slice(1) << 0 << 0 << 0        << endr
                  << 0 << 0 << 0        << endr
                  << 0 << 0 << 13662e-6 << endr;

    I2__.slice(2) << 0  << 0  << 0 << endr
                  << 0  << 0  << 0 << endr
                  << 0  << 0  << 0 << endr;

    I2__.slice(3) << 0  << 0  << 0 << endr
                  << 0  << 0  << 0 << endr
                  << 0  << 0  << 0 << endr;

    Serial PRRP = Serial(4, l2_, lg2_ , m2_, I2__ , {9.8,0, 0}, &fDH_PRRP);


    Serial **R_ = new Serial* [3];
    R_[0] = &_6R1;
    R_[1] = &_6R2;
    R_[2] = &PRRP;

    double b = 0.260;
    double h = 0.008;
    double l = 0.080;
    double d = 0.0039;

    //Matrizes que descrevem a arquitetura do mecanismo
    mat D_ = join_vert(join_vert((mat)eye(3,3),3), zeros(1,3) );
    mat E_ = join_vert( join_diag(join_diag( Rotx(PI/2)*Roty(PI/2), Rotx( PI/2)*Roty(-PI/2)), Roty(-PI/2)*Rotx(PI)) , zeros(1,9) );
    mat F_ = zeros(10,16); F_(9, 13) = 1; F_(9, 14) = 1;
    vec f_ = {0, b-l, h, 0, -(b-l), h, 0, 0, d, -PI/2};
    mat P_ = zeros(6,3);
    mat Q_ = join_horiz(eye(6,6), zeros(6,3) );
    mat S_ = zeros(6,16);

    Parallel Robot = Parallel(3, &P, R_, 3, {3,9,15}, D_, E_, F_, f_, P_, Q_, S_);
    Reference RefObj = Reference(0.12, {0.0, 0.0, 0.480}, {-0.25, -0.25, 0.500});

    /*
    double lambda = 50.0;
    FLControlLaw FL = FLControlLaw(19, lambda*lambda, 2*lambda, &RefObj, &Robot);
    //FLControlLaw FL = FLControlLaw(19, 400.0, 40.0, &r_, &dr_, &d2r_, &Robot);
    Acceleration AC = Acceleration(19, &Robot, &FL);


    vec q0_ = {0.0, 0.0, 0.480,
               9.6321e-01, 0, 1.7529, 0, -1.1452, 0,
               9.6321e-01, 0, 1.7529, 0, -1.1452, 0,
               0, PI/2, 0, 9.6100e-02};
    vec x0_ = join_vert(q0_, zeros(19,1));
    */

    Acceleration AC = Acceleration(16, &Robot, &RefObj);
    vec q0_ = {9.6321e-01, 0, 1.7529, 0, -1.1452, 0,
               9.6321e-01, 0, 1.7529, 0, -1.1452, 0,
               0, PI/2, 0, 9.6100e-02};
    vec x0_ = join_vert(q0_, zeros(16,1));
    

    //Robot.Doit(q0_, zeros(19,1));
    //cout << Robot._q_ << endl;
    //vec aux1_ = {0, b-l, h};
    //vec aux2_ = {0, -(b-l), h};
    //vec aux3_ = {0, 0, d};
    //cout << Rotx( PI/2)*Roty( PI/2)*(*Robot.o__[0]) + aux1_ << endl;
    //cout << Rotx( PI/2)*Roty(-PI/2)*(*Robot.o__[1]) + aux2_ << endl;
    //cout << Roty(-PI/2)*Rotx( PI  )*(*Robot.o__[2]) + aux3_ << endl;

    RK rk = RK("RK8", &AC);
    rk.Doit(0.001, 2*0.12, x0_);
    //for(uint i = 0; i< rk.t_.n_rows; i++) cout << rk.t_(i) << "; " << rk.u__(0,0,i)  << "; " << rk.u__(1,0,i) << "; " << rk.u__(2,0,i) << "; " << endl;
    //for(uint i = 0; i< rk.t_.n_rows; i++) cout << rk.t_(i) << "; " << r_(rk.t_(i))(0) - rk.y__(0,0,i)  << "; " << r_(rk.t_(i))(1) - rk.y__(1,0,i) << "; " << r_(rk.t_(i))(2) - rk.y__(2,0,i) << "; " << endl;
    //for(uint i = 0; i< rk.t_.n_rows; i++) cout << rk.t_(i) << "; " << rk.y__(0,0,i)  << "; " << rk.y__(6,0,i) << "; " << rk.y__(12,0,i) << "; " << endl;

    for(uint i = 0; i< rk.t_.n_rows; i++){
      RefObj.Doit(rk.t_(i));
      //cout << rk.t_(i) << "; " << RefObj.r_(0) - rk.y__(0,0,i)  << "; " << RefObj.r_(1) - rk.y__(1,0,i) << "; " << RefObj.r_(2) - rk.y__(2,0,i) << "; " << endl;
      vec qo0_ = rk.y__(span(0,15) , span(0,0) , span(i,i) );
      vec qo1_ = rk.y__(span(16,31) ,span(0,0) , span(i,i) );
      Robot.Doit(join_vert(RefObj.r_,  qo0_), join_vert(RefObj.dr_, qo1_ ) );
      cout << rk.t_(i) << "; " << norm(Robot._q_) << "; " << 0 << "; " << 0 << "; " <<  endl;
    }
    return 0;
}