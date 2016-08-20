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
    double lg1 = 0.5*0.12;
    double lg2 = 0.5*0.15;
    double m1 = 0.143;
    double m2 = 0.171;
    double Jz1 = 171.6e-6;
    double Jz2 = 320.6e-6;

    double sigma1 = 0.1;
    double sigma2 = 0.1;
    double sigma3 = 0.1;
    double sigma4 = 0.1;
    double sigma5 = 0.3;
    double sigma6 = 0.3;

    cube I__; I__.zeros(3,3,2);
    I__.slice(0) << 0 << 0   << 0   << endr
                 << 0 << Jz1 << 0   << endr
                 << 0 << 0   << Jz1 << endr;

    I__.slice(1) << 0 << 0   << 0   << endr
                 << 0 << Jz2 << 0   << endr
                 << 0 << 0   << Jz2 << endr;

    Mecanismo P = Mecanismo(2);
    Serial RR1 = Serial(2, {0.12, 0.15}, {lg1, lg2},{m1, m2}, I__ , {0, -9.8,0}, &fDH_RR);
    Serial RR2 = Serial(2, {0.12, 0.15}, {lg1, lg2},{m1, m2}, I__ , {0, -9.8,0}, &fDH_RR);
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

    uint nPI = 6;
    uint nCases = (uint)pow(2,nPI);

    Mecanismo **P_ = new Mecanismo* [nCases];
    Serial **RR1_ = new Serial* [nCases];
    Serial **RR2_ = new Serial* [nCases];
    Serial ***RR__ = new Serial** [nCases];
    Parallel **Robot_ = new Parallel* [nCases];

    uint cont = 0;
    for(int i1=-1; i1<2; i1 = i1 + 2){
        for(int i2=-1; i2<2; i2 = i2 + 2){
            for(int i3=-1; i3<2; i3 = i3 + 2){
                for(int i4=-1; i4<2; i4 = i4 + 2){
                    for(int i5=-1; i5<2; i5 = i5 + 2){
                        for(int i6=-1; i6<2; i6 = i6 + 2){
                            I__(1,1,0) = Jz1*(1+sigma5*i5);
                            I__(2,2,0) = Jz1*(1+sigma5*i5);
                            I__(1,1,1) = Jz2*(1+sigma6*i6);
                            I__(2,2,1) = Jz2*(1+sigma6*i6);
                            P_[cont] = new Mecanismo(2);
                            RR1_[cont] = new Serial(2, {0.12, 0.15}, {lg1*(1+sigma1*i1), lg2*(1+sigma2*i2)}, {m1*(1+sigma3*i3), m2*(1+sigma4*i4)}, I__ , {0, -9.8,0}, &fDH_RR);
                            RR2_[cont] = new Serial(2, {0.12, 0.15}, {lg1*(1+sigma1*i1), lg2*(1+sigma2*i2)}, {m1*(1+sigma3*i3), m2*(1+sigma4*i4)}, I__ , {0, -9.8,0}, &fDH_RR);
                            RR__[cont] = new Serial* [2];
                            RR__[cont][0] = RR1_[cont];
                            RR__[cont][1] = RR2_[cont];
                            Robot_[cont] = new Parallel(2, P_[cont], RR__[cont], 2, {2,4}, D_, E_, F_, f_);
                            cont++;
                        }
                    }
                }
            }
        }
    }

    //Plotar área de trabalho
    double lx = 0.23;
    double ly = 0.27;
    double dl = 0.001;
    uint nx = (lx/dl);
    uint ny = (ly/dl);
    Mat<int> M; M.zeros(ny,nx);
    //mat M; M.zeros(ny,nx);
    field<mat> fMh_(ny,nx);
    field<vec> fgh_(ny,nx);
    field<vec> fa1_(ny,nx);
    field<vec> fa2_(ny,nx);
    field<vec> fa12_(ny,nx);
    for(uint i=0; i<ny; i++){
        for(uint j=0; j<nx; j++){
            fMh_(i,j).zeros(2,2);
            fgh_(i,j).zeros(2);
            fa1_(i,j).zeros(2);
            fa2_(i,j).zeros(2);
            fa12_(i,j).zeros(2);
        }
    }
    double r = 0.07;
    double x0 = 0.0;
    double y0 = 0.17;

    uint rows = ny;
    uint cols = nx;

    vec q0_ = {0.823167, 1.81774, 0.823167, 1.81774};

    GNR2 gnr2 = GNR2("RK6", &Robot, 1e-6, 30);

    double x;
    double y;
    mat A2_;
    mat Delta_; Delta_.zeros(2,2);
    vec delta_; delta_.zeros(2);
    vec delta1_; delta1_.zeros(2);
    vec delta2_; delta2_.zeros(2);
    vec delta12_; delta12_.zeros(2);
    vec v1_ = {1,0};
    vec v2_ = {0,1};
    vec v12_ = {1,1};
    vec a1_;
    vec a2_;
    vec a12_;
    for(uint i=0; i<rows; i++){
        for(uint j=0; j<cols; j++){
            x = j*dl+0.5*dl;
            y = i*dl+0.5*dl;
            gnr2.Doit(q0_, {x, y});
            if(gnr2.convergiu){
                A2_ = join_horiz(Robot.Ah_, join_horiz(Robot.Ao_.col(1), Robot.Ao_.col(3)) );
                if(abs(det(Robot.Ao_)) < 1.6*1e-4 || abs(det(A2_)) < 1e-3 ) M(i,j) = 2;
                else{
                    M(i,j) = 1;
                    Robot.Doit(Robot.q0_, Robot.C_*v1_);
                    fMh_(i,j) = Robot.dy->Mh_;
                    fgh_(i,j) = Robot.dy->gh_;
                    fa1_(i,j) = Robot.dy->vh_;
                    Robot.Doit(Robot.q0_, Robot.C_*v2_);
                    fa2_(i,j) = Robot.dy->vh_;
                    Robot.Doit(Robot.q0_, Robot.C_*v12_);
                    fa12_(i,j) = Robot.dy->vh_ - fa1_(i,j) - fa2_(i,j);
                    for (uint k = 0; k < nCases; k++){
                        Robot_[k]->Doit(Robot.q0_, Robot.C_*v1_);
                        Delta_ = arma::max(Delta_, abs(solve(Robot_[k]->dy->Mh_, fMh_(i,j)) - eye(2,2)));
                        delta_ = arma::max(delta_, abs(solve(Robot_[k]->dy->Mh_, fgh_(i,j) - Robot_[k]->dy->gh_)));
                        a1_ = Robot_[k]->dy->vh_;
                        delta1_ = arma::max(delta1_, abs(solve(Robot_[k]->dy->Mh_, fa1_(i,j) - a1_)));
                        Robot_[k]->Doit(Robot.q0_, Robot.C_*v2_);
                        a2_ = Robot_[k]->dy->vh_;
                        delta2_ = arma::max(delta2_, abs(solve(Robot_[k]->dy->Mh_, fa2_(i,j) - a2_)));
                        Robot_[k]->Doit(Robot.q0_, Robot.C_*v12_);
                        a12_ = Robot_[k]->dy->vh_ - a1_ - a2_;
                        delta12_ = arma::max(delta12_, abs(solve(Robot_[k]->dy->Mh_, fa12_(i,j) - a12_)));
                    }
                }
                //M(i,j) = max(log(1 + 1.0/abs(det(Robot.Ao_))), log(1 + 1.0/abs(det(A2_))));
            }
            gnr2.convergiu = false;
            if( ((x - x0)*(x - x0) + (y - y0)*(y - y0) <= r*r) && (M(i,j) != 2) )
                M(i,j) = 3;
        }
    }

    //for(uint i=0; i<rows; i++){
    //    for(uint j=0; j<cols; j++){
    //        cout << M(i,j) << ";" ;
    //        if(j==cols-1) cout << endl;
    //    }
    //}

    cout << Delta_ << endl;
    cout << delta_ << endl;
    cout << delta1_ << endl;
    cout << delta2_ << endl;
    cout << delta12_ << endl;


    return 0;
}