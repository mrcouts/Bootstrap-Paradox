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
#include "GNR2.h"

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
    double lx = 0.23;
    double ly = 0.27;
    //double dl = 0.00025;
    double dl = 0.0005;
    uint nx = (lx/dl);
    uint ny = (ly/dl);
    Mat<int> M; M.zeros(ny,nx);
    double r = 0.07;
    double x0 = 0.0;
    double y0 = 0.17;

    uint rows = ny;
    uint cols = nx;

    vec q0_ = {0.823167, 1.81774, 0.823167, 1.81774};

    GNR2 gnr2 = GNR2("RK6", &Robot, 1e-6, 30);
    //gnr2.Doit(q0_, {0.05,0.08});
    //cout << gnr2.convergiu << endl;
    //cout << gnr2.x_ << endl;
    //cout << gnr2.res_ << endl;
    //cout << gnr2.n << endl;

    double x;
    double y;
    mat A2_;
    for(uint i=0; i<rows; i++){
        for(uint j=0; j<cols; j++){
            x = j*dl+0.5*dl;
            y = i*dl+0.5*dl;
            gnr2.Doit(q0_, {x, y});
            if(gnr2.convergiu){
                q0_ = gnr2.x_;
                A2_ = join_horiz(Robot.Ah_, join_horiz(Robot.Ao_.col(1), Robot.Ao_.col(3)) );
                if(abs(det(Robot.Ao_)) < 1.6*1e-4 || abs(det(A2_)) < 1e-10 ) M(i,j) = 2;
                else M(i,j) = 1;
            }
            gnr2.convergiu = false;
            if( ((x - x0)*(x - x0) + (y - y0)*(y - y0) <= (r+dl)*(r+dl) ) && ((x - x0)*(x - x0) + (y - y0)*(y - y0) >= (r-dl)*(r-dl) ) && (M(i,j) != 2) )
                M(i,j) = 3;
        }
    }

    for(uint i=0; i<rows; i++){
        for(uint j=0; j<cols; j++){
            cout << M(i,j) << ";" ;
            if(j==cols-1) cout << endl;
        }
    }


    return 0;
}