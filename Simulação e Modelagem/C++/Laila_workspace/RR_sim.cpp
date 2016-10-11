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
#include "GNR2.h"

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

    //Plotar área de trabalho
    double lx = 0.42;
    double ly = 0.42;
    //double dl = 0.00025;
    double dl = 0.00125;
    uint nx = (lx/dl);
    uint ny = (ly/dl);
    Mat<int> M; M.zeros(ny,nx);
    double r = 0.37;
    double x0 = 0.0;
    double y0 = 0.0;

    uint rows = ny;
    uint cols = nx;

    vec q0_ = {9.6321e-01, 0, 1.7529, 0, -1.1452, 0,
               9.6321e-01, 0, 1.7529, 0, -1.1452, 0,
               0, PI/2, 0, 9.6100e-02};

    GNR2 gnr2 = GNR2("RK6", &Robot, 1e-6, 30);
    //gnr2.Doit(q0_, {0.05,0.08});
    //cout << gnr2.convergiu << endl;
    //cout << gnr2.x_ << endl;
    //cout << gnr2.res_ << endl;
    //cout << gnr2.n << endl;

    double x;
    double y;
    //mat A2_;
    for(uint i=0; i<rows; i++){
        for(uint j=0; j<cols; j++){
            x = j*dl+0.5*dl;
            y = i*dl+0.5*dl;
            gnr2.Doit(q0_, {x, y, 0.480});
            if(gnr2.convergiu){
                q0_ = gnr2.x_;
                //A2_ = join_horiz(Robot.Ah_, join_horiz(Robot.Ao_.col(1), Robot.Ao_.col(3)) );
                //if(abs(det(Robot.Ao_)) < 1.6*1e-4 || abs(det(A2_)) < 1e-10 ) M(i,j) = 2;
                //if(abs(det(Robot.Ao_)) < 0.2*1e-5) M(i,j) = 2;
                if(0) M(i,j) = 2;
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