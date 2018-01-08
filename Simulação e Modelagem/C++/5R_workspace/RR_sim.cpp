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

    cube I1__; I1__.zeros(3,3,2);
    I1__.slice(0) << 0 << 0                       << 0          << endr
                  << 0 << 107.307e-6 + 146.869e-6 << 0          << endr
                  << 0 << 0                       << 107.307e-6 + 146.869e-6 << endr;

    I1__.slice(1) << 0 << 0        << 0        << endr
                  << 0 << 438.0e-6 << 0        << endr
                  << 0 << 0        << 438.0e-6 << endr;

    cube I2__; I2__.zeros(3,3,2);
    I2__.slice(0) << 0 << 0                        << 0          << endr
                  << 0 << 107.307e-6 + 188.738e-6  << 0          << endr
                  << 0 << 0                        << 107.307e-6 + 188.738e-6  << endr;

    I2__.slice(1) << 0 << 0          << 0          << endr
                  << 0 << 301.679e-6 << 0          << endr
                  << 0 << 0          << 301.679e-6 << endr;

    Serial RR1 = Serial(2, {0.12, 0.16}, {0.06, 0.078},{0.062, 0.124}, I1__ , {0, 0, 9.8}, &fDH_RR);
    Serial RR2 = Serial(2, {0.12, 0.16}, {0.06, 0.058},{0.062, 0.097}, I2__ , {0, 0, 9.8}, &fDH_RR);
    Serial **RR_ = new Serial* [2];
    RR_[0] = &RR1;
    RR_[1] = &RR2;

    //Matrizes que descrevem a arquitetura do mecanismo
    double l0 = 0.05;
    mat D_ = join_vert((mat)eye(2,2),2);
    mat E_ = join_diag( Roty(0)(span(0,1),span(0,2)), Roty(PI)(span(0,1),span(0,2)) );
    mat F_ = zeros(4,4);
    vec f_ = {l0,0,-l0,0};

    Parallel Robot = Parallel(2, &P, RR_, 2, {2,4}, D_, E_, F_, f_);
    Reference RefObj = Reference(0.12, {0.08, 0.16}, {-0.08, 0.4});

    //Plotar área de trabalho
    uint nx = 96.0;
    uint ny = 56.0;
    double lx = 0.24;
    double ly = 0.28;
    double xi = -lx;
    double xf = lx;
    double yi = 0.0;
    double yf = ly;
    double dx = (xf-xi)/(nx-1);
    double dy = (yf-yi)/(ny-1);
    double dl = 0.5*(dx+dy);
    
    Mat<int> M; M.zeros(nx,ny);
    field<mat> fZ_(nx,ny);
    field<mat> fMh_(nx,ny);
    field<vec> fgh_(nx,ny);
    field<vec> fa1_(nx,ny);
    field<vec> fa2_(nx,ny);
    field<vec> fa12_(nx,ny);
    for(uint i=0; i<nx; i++){
        for(uint j=0; j<ny; j++){
            fZ_(i,j).zeros(2,2);
            fMh_(i,j).zeros(2,2);
            fgh_(i,j).zeros(2);
            fa1_(i,j).zeros(2);
            fa2_(i,j).zeros(2);
            fa12_(i,j).zeros(2);
        }
    }
    vec v1_ = {1,0};
    vec v2_ = {0,1};
    vec v12_ = {1,1};
    double r = 0.07;
    double x0 = 0.0;
    double y0 = 0.17;

    uint rows = nx;
    uint cols = ny;

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
            x = xi + i*dx;
            y = yi + j*dy;
            gnr2.Doit(q0_, {x, y});
            if(gnr2.convergiu){
                q0_ = gnr2.x_;
                A2_ = join_horiz(Robot.Ah_, join_horiz(Robot.Ao_.col(1), Robot.Ao_.col(3)) );
                if(abs(det(Robot.Ao_)) < 1.6*1e-6 || abs(det(A2_)) < 1e-11 ) M(i,j) = 2;
                else{ 
                	M(i,j) = 1;
                	Robot.Doit(Robot.q0_, join_vert(v1_, -solve(Robot.Ao_, Robot.Ah_*v1_)) );
                	fZ_(i,j) = Robot.Z_;
                    fMh_(i,j) = Robot.dy->Mh_;
                    fgh_(i,j) = Robot.dy->gh_;
                    fa1_(i,j) = Robot.dy->vh_;
                    Robot.Doit(Robot.q0_, Robot.C_*v2_);
                    fa2_(i,j) = Robot.dy->vh_;
                    Robot.Doit(Robot.q0_, Robot.C_*v12_);
                    fa12_(i,j) = Robot.dy->vh_ - fa1_(i,j) - fa2_(i,j);
                }
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


    fZ_.save("fZ_field");
    fMh_.save("fMh_field");
    fgh_.save("fgh_field");
    fa1_.save("fa1_field");
    fa2_.save("fa2_field");
    fa12_.save("fa12_field");


    return 0;
}