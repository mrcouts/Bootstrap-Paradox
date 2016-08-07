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
    double r = 0.08;
    double x0 = 0.0;
    double y0 = 0.16;
    double w = 1/r;
    return {x0+r*cos(w*t), y0+r*sin(w*t)};
}

vec dr_(double t){
    double r = 0.08;
    double w = 1/r;
    return {-w*r*sin(w*t), w*r*cos(w*t)};
}

vec d2r_(double t){
    double r = 0.08;
    double w = 1/r;
    return {-w*w*r*cos(w*t), -w*w*r*sin(w*t)};
}

Dy* (dy_comp)(vec q0_, vec q1_){
    Dy *dy;
    dy = new Dy(2);
    return dy;
}

int main(){

    Mecanismo P = Mecanismo(2);

        cube I__; I__.zeros(3,3,2);
    I__.slice(0) << 0 << 0        << 0        << endr
                 << 0 << 171.6e-6 << 0        << endr
                 << 0 << 0        << 171.6e-6 << endr;

    I__.slice(1) << 0 << 0        << 0        << endr
                 << 0 << 320.6e-6 << 0        << endr
                 << 0 << 0        << 320.6e-6 << endr;

    Serial RR1 = Serial(2, {0.12, 0.15}, {0.5*0.12, 0.5*0.15},{0.143, 0.171}, I__ , {0, -9.8,0}, &fDH_RR);
    Serial RR2 = Serial(2, {0.12, 0.15}, {0.5*0.12, 0.5*0.15},{0.143, 0.171}, I__ , {0, -9.8,0}, &fDH_RR);
    Serial **RR_ = new Serial* [2];
    RR_[0] = &RR1;
    RR_[1] = &RR2;

    //Matrizes que descrevem a arquitetura do mecanismo
    mat D_ = join_vert((mat)eye(2,2),2);
    mat E_ = join_diag( Roty(0)(span(0,1),span(0,2)), Roty(PI)(span(0,1),span(0,2)) );
    mat F_ = zeros(4,4);
    vec f_ = {0.05,0,-0.05,0};

    Parallel Robot = Parallel(2, &P, RR_, 2, {2,4}, D_, E_, F_, f_);
    Reference RefObj = Reference(0.12, {0.08, 0.16}, {-0.08, 0.4});

    //Criando lista de Objetos tipo trajetória
    Reference **RefObjVec = new Reference* [5];
    for(uint i=0; i<5; i++) RefObjVec[i] = new Reference(0.12, {0.08, 0.16}, {-0.08, 0.13+0.1*i});
    RefObjVec[2]->Doit(0.05);
    cout << RefObjVec[2]->d2r_ << endl;
    delete[] RefObjVec; 


    //Simulação dinâmica
    /*
    double lambda = 50.0;
    FLControlLaw FL = FLControlLaw(19, lambda*lambda, 2*lambda, &RefObj, &Robot);
    //FLControlLaw FL = FLControlLaw(6, lambda*lambda, 2*lambda, &r_, &dr_, &d2r_, &Robot);
    Acceleration AC = Acceleration(6, &Robot, &FL);

    vec q0_ = {0.08,  0.16, 0.305030291698133, 1.86386236511897, 1.45111035931733, 1.41460649673445};
    vec x0_ = join_vert(q0_, zeros(6));
    */

    //Simulação cinemática
    Acceleration AC = Acceleration(4, &Robot, &RefObj);
    vec q0_ = {0.305030291698133, 1.86386236511897, 1.45111035931733, 1.41460649673445};
    vec x0_ = join_vert(q0_, zeros(4));

    RK rk = RK("RK8", &AC);
    //rk.Doit(0.001, 2*0.12, x0_);
    //for(uint i = 0; i< rk.t_.n_rows; i++) cout << rk.t_(i) << "; " << rk.u__(0,0,i)  << "; " << rk.u__(1,0,i) << "; " << endl;
    //for(uint i = 0; i< rk.t_.n_rows; i++) cout << rk.t_(i) << "; " << r_(rk.t_(i))(0) - rk.y__(0,0,i)  << "; " << r_(rk.t_(i))(1) - rk.y__(1,0,i) << "; " << endl;
    //for(uint i = 0; i< rk.t_.n_rows; i++){
    //  RefObj.Doit(rk.t_(i));
    //  //cout << rk.t_(i) << "; " << RefObj.r_(0) - rk.y__(0,0,i)  << "; " << RefObj.r_(1) - rk.y__(1,0,i) << "; " << RefObj.r_(2) - rk.y__(2,0,i) << "; " << endl;
    //  vec qo0_ = rk.y__(span(0,3), span(0,0), span(i,i) );
    //  vec qo1_ = rk.y__(span(4,7), span(0,0), span(i,i) );
    //  Robot.Doit(join_vert(RefObj.r_,  qo0_), join_vert(RefObj.dr_, qo1_ ) );
    //  cout << rk.t_(i) << "; " << norm(Robot._q_, "inf") << "; " << arma::rank(Robot.A_) << "; " <<  endl;
    //}
    return 0;
}