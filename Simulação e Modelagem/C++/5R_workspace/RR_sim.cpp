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

    //Plotar área de trabalho
    double lx = 0.25;
    double ly = 0.28;
    double dl = 0.005;
    uint nx = (lx/dl);
    uint ny = (ly/dl);
    Mat<int> M; M.zeros(ny,nx);
    double r = 0.085;
    double x0 = 0.0;
    double y0 = 0.16;

    uint rows = ny;
    uint cols = nx;

    vec q0_ = {0.823167, 1.81774, 0.823167, 1.81774};
    vec x0_ = join_vert(q0_, zeros(4));

    Reference ***RefObjMat = new Reference** [rows];
    Acceleration ***ACMat = new Acceleration** [rows];
    RK ***rkMat = new RK** [rows];
    for(uint i=0; i<rows; i++){
        RefObjMat[i] = new Reference* [cols];
        ACMat[i] = new Acceleration* [cols];
        rkMat[i] = new RK* [cols];
        for(uint j=0; j<cols; j++){
            RefObjMat[i][j] = new Reference(0.020, {0.00, 0.016}, {j*dl+0.5*dl, i*dl+0.5*dl});
            ACMat[i][j] = new Acceleration(4, &Robot, RefObjMat[i][j]);
            rkMat[i][j] = new RK("RK8", ACMat[i][j]);
            rkMat[i][j]->Doit(0.001, 0.024, x0_);
            if( arma::rank(Robot.A_) != 0 && abs(norm(Robot._q_, "inf") ) < 1e-4 ) M(i,j) = 1;
        }
    }

    for(uint i=0; i<rows; i++){
        for(uint j=0; j<cols; j++){
            cout << M(i,j) << ";" ;
            if(j==cols-1) cout << endl;
        }
    }

    for(uint i = 0; i < rows; i++){
        delete [] RefObjMat[i];
        delete [] ACMat[i];
        delete [] rkMat[i];
    }
    delete [] RefObjMat;
    delete [] ACMat;
    delete [] rkMat;

    return 0;
}