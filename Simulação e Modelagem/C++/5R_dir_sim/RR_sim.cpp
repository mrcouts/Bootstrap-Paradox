#include <iostream>
#include <iomanip>
#include <string>
#include "SomeUtilities.h"
#include "Serial.h"
#include "RR.h"
#include "FLControlLaw.h"
#include "Acceleration.h"
#include "RK.h"
#include "Parallel.h"

vec r_(double t){
    double r = 0.05;
    double x0 = 0.0;
    double y0 = 0.158;
    double w = 2*PI;
    if(t<8.0){
        return {x0-r*sin(w*t), y0-r*cos(w*t)};}
    else
        return {x0-r*sin(w*0), y0-r*cos(w*0)};
}

vec dr_(double t){
    double r = 0.05;
    double w = 2*PI;
    if(t<8.0){
        return {-w*r*cos(w*t), +w*r*sin(w*t)};}
    else
        return {0,0};
}

vec d2r_(double t){
    double r = 0.05;
    double w = 2*PI;
    if(t<8.0){
        return {+w*w*r*sin(w*t), w*w*r*cos(w*t)};}
    else
        return {0,0};
}

Dy* (dy_comp)(vec q0_, vec q1_){
    Dy *dy;
    dy = new Dy(2);
    return dy;
}

int main(){

    Mecanismo P = Mecanismo(2);
    Mecanismo _P= Mecanismo(2);

//    double l1  = 0.12;
//    double l2  = 0.16;
//    double lg1 = 0.06;
//    double lg2 = 0.68;
//    double m1 = 0.062;
//    double m2 = 0.110;
//    double Jz1 = 200.6e-6;
//    double Jz2 = 750.6e-6;
//
//    double l1  = 0.12;
//    double l2  = 0.16;
//    double lg1 = 0.06;
//    double lg2 = 0.68;
//    double m1 = 0.062;
//    double m2 = 0.110;
//    double Jz1 = 200.6e-6;
//    double Jz2 = 750.6e-6;
//
//    arma_rng::set_seed_random();
//    //vec sigma_ = ((vec){0.3, 0.1, 0.2, 0.1, 0.5, 0.1}) % sign(randn(6))*1.8;
//    //vec sigma_ = sign(randn(6))*0.5;
//    //vec sigma_ = {-0.15,0.15,0.15,0.15,0.15,0.15};
    vec sigma_ = {-1.0,1.0,1.0,1.0,1.0,1.0};
    vec coef_  = ones(6) + 0.2*sigma_;
//
//    cube I__; I__.zeros(3,3,2);
//    I__.slice(0) << 0 << 0   << 0   << endr
//                 << 0 << Jz1 << 0   << endr
//                 << 0 << 0   << Jz1 << endr;
//
//    I__.slice(1) << 0 << 0   << 0   << endr
//                 << 0 << Jz2 << 0   << endr
//                 << 0 << 0   << Jz2 << endr;
//
//    cube _I__; _I__.zeros(3,3,2);
//    _I__.slice(0) << 0 << 0           << 0            << endr
//                 << 0 << coef_(4)*Jz1 << 0            << endr
//                 << 0 << 0            << coef_(4)*Jz1 << endr;
//
//    _I__.slice(1) << 0 << 0            << 0            << endr
//                  << 0 << coef_(5)*Jz2 << 0            << endr
//                  << 0 << 0            << coef_(5)*Jz2 << endr;
//

    cube I1__; I1__.zeros(3,3,2);
    I1__.slice(0) << 0 << 0                       << 0          << endr
                  << 0 << 107.307e-6 + 146.869e-6 << 0          << endr
                  << 0 << 0                       << 107.307e-6 + 146.869e-6 << endr;

    I1__.slice(1) << 0 << 0        << 0        << endr
                  << 0 << 576.9e-6 << 0        << endr
                  << 0 << 0        << 576.9e-6 << endr;

    cube I2__; I2__.zeros(3,3,2);
    I2__.slice(0) << 0 << 0                        << 0          << endr
                  << 0 << 107.307e-6 + 188.738e-6  << 0          << endr
                  << 0 << 0                        << 107.307e-6 + 188.738e-6  << endr;

    I2__.slice(1) << 0 << 0          << 0          << endr
                  << 0 << 301.679e-6 << 0          << endr
                  << 0 << 0          << 301.679e-6 << endr;

    Serial RR1 = Serial(2, {0.12, 0.16}, {0.06, 0.078},{0.062, 0.124}, {1.226e-4, 0}, {0*4.357e-2, 0}, I1__ , {0, 0, 9.8}, &fDH_RR);
    Serial RR2 = Serial(2, {0.12, 0.16}, {0.06, 0.058},{0.062, 0.097}, {2.039e-4, 0}, {0*4.992e-2, 0}, I2__ , {0, 0, 9.8}, &fDH_RR);
    Serial **RR_ = new Serial* [2];
    RR_[0] = &RR1;
    RR_[1] = &RR2;

    //Serial _RR1 = Serial(2, {l1, l2}, {coef_(0)*lg1, coef_(1)*lg2},{coef_(2)*m1, coef_(3)*m2}, _I__ , {0, -0*9.8,0}, &fDH_RR);
    //Serial _RR2 = Serial(2, {l1, l2}, {coef_(0)*lg1, coef_(1)*lg2},{coef_(2)*m1, coef_(3)*m2}, _I__ , {0, -0*9.8,0}, &fDH_RR);
    Serial _RR1 = Serial(2, {0.12, 0.16}, {coef_(0)*0.06, coef_(1)*0.078},{coef_(2)*0.062, coef_(3)*0.124}, {coef_(0)*1.226e-4, 0}, {4.357e-2, 0}, coef_(4)*I1__ , {0, 0, 9.8}, &fDH_RR);
    Serial _RR2 = Serial(2, {0.12, 0.16}, {coef_(0)*0.06, coef_(1)*0.058},{coef_(2)*0.062, coef_(3)*0.097}, {coef_(0)*2.039e-4, 0}, {4.992e-2, 0}, coef_(4)*I2__ , {0, 0, 9.8}, &fDH_RR);
    Serial **_RR_ = new Serial* [2];
    _RR_[0] = &_RR1;
    _RR_[1] = &_RR2;

    //Matrizes que descrevem a arquitetura do mecanismo
    mat D_ = join_vert((mat)eye(2,2),2);
    mat E_ = join_diag( Roty(0)(span(0,1),span(0,2)), Roty(PI)(span(0,1),span(0,2)) );
    mat F_ = zeros(4,4);
    vec f_ = {0.05,0,-0.05,0};

    Parallel Robot = Parallel(2, &P, RR_, 2, {2,4}, D_, E_, F_, f_);
    Parallel _Robot= Parallel(2,&_P,_RR_, 2, {2,4}, D_, E_, F_, f_);
    Reference RefObj = Reference(1.2, {0.07, 0.17}, {0.07, 0.17});

    //Simulação dinâmica
    double lambda = 60.0;
    
    
    //double eta = 93.5768;
    //mat Lambda_; Lambda_.zeros(2,2);
    //Lambda_ << 2.6206e+02 << 1.9077e+02 << endr
    //        << 0.0        << 1.7379e+02 << endr;
    //vec gamma_ = {6.8903, 7.4243};

    //vec eta_ = {19.7503, 28.7944};
    //cube Lambda__; Lambda__.zeros(2,2,2);
    //Lambda__.slice(0) << 8.4176e+01 << 8.3400e+01 << endr
    //                  << 0          << 5.4837e+01 << endr;
    //Lambda__.slice(1) << 1.0149e+02 << 9.5498e+01 << endr
    //                  << 0          << 6.8895e+01 << endr;
    //mat Gamma_; Gamma_.zeros(2,2);
    //Gamma_ << 2.1547 << 1.5500 << endr
    //       << 2.7074 << 2.2725 << endr;
    vec eta_ = {7.5314e+01, 1.0840e+02};
    cube Lambda__; Lambda__.zeros(2,2,2);
    Lambda__.slice(0) << 37.0403 << 22.7291 << endr
                      << 0       << 17.5189 << endr;
    Lambda__.slice(1) << 54.0692 << 35.5080 << endr
                      << 0       << 28.4984 << endr;
    mat Gamma_; Gamma_.zeros(2,2);
    Gamma_ << 1.2525 << 0.9153 << endr
           << 1.6833 << 1.6495 << endr;

    /*
    //SMCLaw SMC = SMCLaw(2, lambda, eta, K_, k_, 100.0, &RefObj, &Robot);
    //SMCLaw SMC = SMCLaw(2, lambda, eta, Lambda_, gamma_, 20.0, &r_, &dr_, &d2r_, &Robot);
    SMCLaw SMC = SMCLaw(2, lambda, eta_, Lambda__, Gamma_, 100.0, &r_, &dr_, &d2r_, &Robot);
    //SMCLaw SMC = SMCLaw(2, lambda, eta_, Lambda__, Gamma_, 20.0, &RefObj, &Robot);
    //SMCLaw SMC = SMCLaw(2, lambda, 20.0, 20.0, &RefObj, &Robot);
    //SMCLaw SMC = SMCLaw(2, lambda, 40.0, 40.0, &r_, &dr_, &d2r_, &Robot);
    Acceleration AC = Acceleration(6, &_Robot, &SMC);
    */

    //FLControlLaw FL = FLControlLaw(2, lambda*lambda, 2*lambda, &RefObj, &Robot);
    FLControlLaw FL = FLControlLaw(2, lambda*lambda, 2*lambda, &r_, &dr_, &d2r_, &Robot);
    Acceleration AC = Acceleration(6, &_Robot, &FL);

   
    
    

    vec q0_ = {0.0,  0.108,  0.5375391526183005, 2.308800210309811,  0.5375391526183005, 2.308800210309811};
    vec dq0_ = {-0.31415927, 0.0, 3.38559762, -1.10585113, -3.38559762, 1.10585113};
    vec x0_ = join_vert(q0_, dq0_);

    

    //Simulação cinemática
    /*
    Acceleration AC = Acceleration(4, &Robot, &RefObj);

    vec q0_ = {0.305030291698133, 1.86386236511897, 1.45111035931733, 1.41460649673445};
    vec x0_ = join_vert(q0_, zeros(4));
    */

    RK rk = RK("RK6", &AC);
    rk.Doit(0.0005, 3.0, x0_);
    double t;
    vec ref_; ref_.zeros(2);
    vec dref_; dref_.zeros(2);
    for(uint i = 0; i< rk.t_.n_rows; i++){
        t = rk.t_(i);
        cout << t << "; " << rk.u__(0,0,i)  << "; " << rk.u__(1,0,i) << "; ";
        ref_ = r_(t);
        dref_ = dr_(t);
        //RefObj.Doit(rk.t_(i));
        //ref_ = RefObj.r_;
        //dref_ = RefObj.dr_;
        cout << ref_(0) - rk.y__(0,0,i)  << "; " << ref_(1) - rk.y__(1,0,i) << "; ";
        cout << -( dref_(0) - rk.y__(6,0,i) + lambda*(ref_(0) - rk.y__(0,0,i)) ) << "; " << -( dref_(1) - rk.y__(7,0,i) + lambda*(ref_(1) - rk.y__(1,0,i))) << "; " << endl;
    }
    return 0;
}   