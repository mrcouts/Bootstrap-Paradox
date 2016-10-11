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
    double lg1 = 0.158;
    double lg4 = 0.175;
    double m1 = 0.100;
    double m4 = 0.060;
    double J1 = 1165e-6;
    double J4 = 672.8e-6;

    double lg22 = 0.247;
    double lg42 = 0.100;
    double m12 = 0.738;
    double m22 = 0.495;
    double m42 = 0.128;
    double J22 = 13662e-6;

    double m0 = 0.050;

    double sigma1 = 0.15;
    double sigma2 = 0.15;
    double sigma3 = 0.15;
    double sigma4 = 0.15;
    double sigma5 = 0.15;
    double sigma6 = 0.15;
    double sigma7 = 0.15;
    double sigma8 = 0.15;
    double sigma9 = 0.15;
    double sigma10 = 0.15;
    double sigma11 = 0.15;
    double sigma12 = 0.15;
    double sigma13 = 0.15;

    vec l_ = {0.305, 0.045, 0.045 , 0.372, 0.0125, 0.000};
    vec lg_= {lg1  , 0.025, 0.0116, lg4  , 0.0042, 0.000};
    vec m_ = {m1   , 0.033, 0.034 , m4   , 0.023 , 0.000};

    cube I__; I__.zeros(3,3,6);
    I__.slice(0) << 0.99*J1 << 0  << 0       << endr
                 << 0       << J1 << 0       << endr
                 << 0       << 0  << 12.7e-6 << endr;

    I__.slice(1) << 8.0e-6 << 0      << 0      << endr
                 << 0      << 1.5e-6 << 0      << endr
                 << 0      << 0      << 8.2e-6 << endr;

    I__.slice(2) << 7.9e-6 << 0      << 0      << endr
                 << 0      << 8.2e-6 << 0      << endr
                 << 0      << 0      << 1.2e-6 << endr;

    I__.slice(3) << J4 << 0      << 0  << endr
                 << 0  << 1.7e-6 << 0  << endr
                 << 0  << 0      << J4 << endr;

    I__.slice(4) << 1.6e-6 << 0      << 0      << endr
                 << 0      << 1.3e-6 << 0      << endr
                 << 0      << 0      << 1.2e-6 << endr;

    I__.slice(5) << 0 << 0 << 0 << endr
                 << 0 << 0 << 0 << endr
                 << 0 << 0 << 0 << endr;

    vec l2_ = {0.000, 0.380, 0.000, 0.200};
    vec lg2_= {0.000, lg22 , 0.000, lg42 };
    vec m2_ = {m12  , m22  , 0.000, m42  };

    cube I2__; I2__.zeros(3,3,4);
    I2__.slice(0) << 0  << 0  << 0 << endr
                  << 0  << 0  << 0 << endr
                  << 0  << 0  << 0 << endr;

    I2__.slice(1) << 0 << 0 << 0   << endr
                  << 0 << 0 << 0   << endr
                  << 0 << 0 << J22 << endr;

    I2__.slice(2) << 0  << 0  << 0 << endr
                  << 0  << 0  << 0 << endr
                  << 0  << 0  << 0 << endr;

    I2__.slice(3) << 0  << 0  << 0 << endr
                  << 0  << 0  << 0 << endr
                  << 0  << 0  << 0 << endr;

    Mecanismo P = Mecanismo(3); P.dy->Mh_ = m0*eye(3,3);
    Serial _6R1 = Serial(6, l_, lg_ , m_, I__ , {0, 9.8,0}, &fDH_6R);
    Serial _6R2 = Serial(6, l_, lg_ , m_, I__ , {0, 9.8,0}, &fDH_6R);
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
    mat Pe_ = zeros(6,3);
    mat Q_ = join_horiz(eye(6,6), zeros(6,3) );
    mat S_ = zeros(6,16);

    Parallel Robot = Parallel(3, &P, R_, 3, {3,9,15}, D_, E_, F_, f_, Pe_, Q_, S_);

    uint nPI = 13;
    uint nCases = (uint)pow(2,nPI);

    Mecanismo **P_ = new Mecanismo* [nCases];
    Serial **_6R1_ = new Serial* [nCases];
    Serial **_6R2_ = new Serial* [nCases];
    Serial **PRRP_ = new Serial* [nCases];
    Serial ***RR__ = new Serial** [nCases];
    Parallel **Robot_ = new Parallel* [nCases];

    uint cont = 0;
    for(int i1=-1; i1<2; i1 = i1 + 2){
        for(int i2=-1; i2<2; i2 = i2 + 2){
            for(int i3=-1; i3<2; i3 = i3 + 2){
                for(int i4=-1; i4<2; i4 = i4 + 2){
                    for(int i5=-1; i5<2; i5 = i5 + 2){
                        for(int i6=-1; i6<2; i6 = i6 + 2){
                            for(int i7=-1; i7<2; i7 = i7 + 2){
                                for(int i8=-1; i8<2; i8 = i8 + 2){
                                    for(int i9=-1; i9<2; i9 = i9 + 2){
                                        for(int i10=-1; i10<2; i10 = i10 + 2){
                                            for(int i11=-1; i11<2; i11 = i11 + 2){
                                                for(int i12=-1; i12<2; i12 = i12 + 2){
                                                    for(int i13=-1; i13<2; i13 = i13 + 2){
                                                        lg_(0) = lg1*(1+sigma1*i1);
                                                        lg_(3) = lg4*(1+sigma2*i2);
                                                        m_(0) = m1*(1+sigma3*i3);
                                                        m_(3) = m4*(1+sigma4*i4);
                                                        I__(0,0,0) = 0.99*J1*(1+sigma5*i5);
                                                        I__(1,1,0) =      J1*(1+sigma5*i5);
                                                        I__(2,2,3) =      J4*(1+sigma6*i6);
                                                        lg2_(1) = lg22*(1+sigma7*i7);
                                                        lg2_(3) = lg42*(1+sigma8*i8);
                                                        m2_(0) = m12*(1+sigma9*i9);
                                                        m2_(1) = m22*(1+sigma10*i10);
                                                        m2_(3) = m42*(1+sigma11*i11);
                                                        I2__(2,2,1) = J22*(1+sigma12*i12);
                                                        P_[cont] = new Mecanismo(3); P_[cont]->dy->Mh_ = m0*(1+sigma13*i13)*eye(3,3);
                                                        _6R1_[cont] = new Serial(6, l_, lg_ , m_, I__ , {0, 9.8,0}, &fDH_6R);
                                                        _6R2_[cont] = new Serial(6, l_, lg_ , m_, I__ , {0, 9.8,0}, &fDH_6R);
                                                        PRRP_[cont] = new Serial(4, l2_, lg2_ , m2_, I2__ , {9.8,0, 0}, &fDH_PRRP);

                                                        RR__[cont] = new Serial* [3];
                                                        RR__[cont][0] = _6R1_[cont];
                                                        RR__[cont][1] = _6R2_[cont];
                                                        RR__[cont][2] = PRRP_[cont];

                                                        Robot_[cont] = new Parallel(3, P_[cont], RR__[cont], 3, {3,9,15}, D_, E_, F_, f_, Pe_, Q_, S_);
                                                        cont++;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //Plotar área de trabalho
    double lx = 0.4;
    double ly = 0.4;
    double dl = 0.2;
    uint nx = (lx/dl);
    uint ny = (ly/dl);
    Mat<int> M; M.zeros(ny,nx);
    //mat M; M.zeros(ny,nx);
    field<mat> fMh_(ny,nx);
    field<vec> fgh_(ny,nx);
    field<vec> fa1_(ny,nx);
    field<vec> fa2_(ny,nx);
    field<vec> fa3_(ny,nx);
    field<vec> fa12_(ny,nx);
    field<vec> fa13_(ny,nx);
    field<vec> fa23_(ny,nx);
    for(uint i=0; i<ny; i++){
        for(uint j=0; j<nx; j++){
            fMh_(i,j).zeros(3,3);
            fgh_(i,j).zeros(3);
            fa1_(i,j).zeros(3);
            fa2_(i,j).zeros(3);
            fa3_(i,j).zeros(3);
            fa12_(i,j).zeros(3);
            fa13_(i,j).zeros(3);
            fa23_(i,j).zeros(3);
        }
    }
    double r = 0.37;
    double x0 = 0.0;
    double y0 = 0.0;

    uint rows = ny;
    uint cols = nx;

    vec q0_ = {9.6321e-01, 0, 1.7529, 0, -1.1452, 0,
               9.6321e-01, 0, 1.7529, 0, -1.1452, 0,
               0, PI/2, 0, 9.6100e-02};

    GNR2 gnr2 = GNR2("RK6", &Robot, 1e-6, 30);

    double x;
    double y;
    //mat A2_;
    mat Delta_; Delta_.zeros(3,3);
    vec delta_; delta_.zeros(3);
    vec delta1_; delta1_.zeros(3);
    vec delta2_; delta2_.zeros(3);
    vec delta3_; delta3_.zeros(3);
    vec delta12_; delta12_.zeros(3);
    vec delta13_; delta13_.zeros(3);
    vec delta23_; delta23_.zeros(3);
    vec v1_ = {1,0,0};
    vec v2_ = {0,1,0};
    vec v3_ = {0,0,1};
    vec v12_ = {1,1,0};
    vec v13_ = {1,0,1};
    vec v23_ = {0,1,1};
    vec a1_ = {0,0,0};
    vec a2_ = {0,0,0};
    vec a3_ = {0,0,0};
    vec a12_ = {0,0,0};
    vec a13_ = {0,0,0};
    vec a23_ = {0,0,0};



    for(uint i=0; i<rows; i++){
        for(uint j=0; j<cols; j++){
            x = j*dl+0.5*dl;
            y = i*dl+0.5*dl;
            gnr2.Doit(q0_, {x, y, 0.480});
            if(gnr2.convergiu){
            	q0_ = gnr2.x_;
                //A2_ = join_horiz(Robot.Ah_, join_horiz(Robot.Ao_.col(1), Robot.Ao_.col(3)) );
                //if(abs(det(Robot.Ao_)) < 1.6*1e-4 || abs(det(A2_)) < 1e-10 ) M(i,j) = 2;
                if( (x - x0)*(x - x0) + (y - y0)*(y - y0) > r*r  ) M(i,j) = 2;
                else{
                    M(i,j) = 1;
                    Robot.Doit(Robot.q0_, join_vert(v1_, -solve(Robot.Ao_, Robot.Ah_*v1_)) );
                    fMh_(i,j) = Robot.dy->Mh_;
                    fgh_(i,j) = Robot.dy->gh_;
                    fa1_(i,j) = Robot.dy->vh_;
                    Robot.Doit(Robot.q0_, Robot.C_*v2_);
                    fa2_(i,j) = Robot.dy->vh_;
                    Robot.Doit(Robot.q0_, Robot.C_*v3_);
                    fa3_(i,j) = Robot.dy->vh_;
                    Robot.Doit(Robot.q0_, Robot.C_*v12_);
                    fa12_(i,j) = Robot.dy->vh_ - fa1_(i,j) - fa2_(i,j);
                    Robot.Doit(Robot.q0_, Robot.C_*v13_);
                    fa13_(i,j) = Robot.dy->vh_ - fa1_(i,j) - fa3_(i,j);
                    Robot.Doit(Robot.q0_, Robot.C_*v23_);
                    fa23_(i,j) = Robot.dy->vh_ - fa2_(i,j) - fa3_(i,j);


                    for (uint k = 0; k < nCases; k++){                      
                        Robot_[k]->Doit(Robot.q0_, Robot.C_*v1_); 
                        Delta_ = arma::max(Delta_, abs(solve(Robot_[k]->dy->Mh_, fMh_(i,j)) - eye(3,3)));
                        delta_ = arma::max(delta_, abs(solve(Robot_[k]->dy->Mh_, fgh_(i,j) - Robot_[k]->dy->gh_)));
                        a1_ = Robot_[k]->dy->vh_;
                        delta1_ = arma::max(delta1_, abs(solve(Robot_[k]->dy->Mh_, fa1_(i,j) - a1_)));



                        Robot_[k]->Doit(Robot.q0_, Robot.C_*v2_);
                        a2_ = Robot_[k]->dy->vh_;
                        delta2_ = arma::max(delta2_, abs(solve(Robot_[k]->dy->Mh_, fa2_(i,j) - a2_)));

                        Robot_[k]->Doit(Robot.q0_, Robot.C_*v3_);
                        a3_ = Robot_[k]->dy->vh_;
                        delta3_ = arma::max(delta3_, abs(solve(Robot_[k]->dy->Mh_, fa3_(i,j) - a3_)));

                        Robot_[k]->Doit(Robot.q0_, Robot.C_*v12_);
                        a12_ = Robot_[k]->dy->vh_ - a1_ - a2_;
                        delta12_ = arma::max(delta12_, abs(solve(Robot_[k]->dy->Mh_, fa12_(i,j) - a12_)));

                        Robot_[k]->Doit(Robot.q0_, Robot.C_*v13_);
                        a13_ = Robot_[k]->dy->vh_ - a1_ - a3_;
                        delta13_ = arma::max(delta13_, abs(solve(Robot_[k]->dy->Mh_, fa13_(i,j) - a13_)));

                        Robot_[k]->Doit(Robot.q0_, Robot.C_*v23_);
                        a23_ = Robot_[k]->dy->vh_ - a2_ - a3_;
                        delta23_ = arma::max(delta23_, abs(solve(Robot_[k]->dy->Mh_, fa23_(i,j) - a23_)));
                    }
                }
                //M(i,j) = max(log(1 + 1.0/abs(det(Robot.Ao_))), log(1 + 1.0/abs(det(A2_))));
            }
            gnr2.convergiu = false;
            //if( ((x - x0)*(x - x0) + (y - y0)*(y - y0) <= r*r) && (M(i,j) != 2) )
            //    M(i,j) = 3;
        }
    }

    //for(uint i=0; i<rows; i++){
    //    for(uint j=0; j<cols; j++){
    //        cout << M(i,j) << ";" ;
    //        if(j==cols-1) cout << endl;
    //    }
    //}

    Delta_.print("Delta_ =");
    delta_.print("delta_ =");
    delta1_.print("delta1_ =");
    delta2_.print("delta2_ =");
    delta3_.print("delta3_ =");
    delta12_.print("delta12_ =");
    delta13_.print("delta13_ =");
    delta23_.print("delta23_ =");

    vec den_ = ones(3) - sum(Delta_,1);
    //cout << den_ << endl;
    den_.print("den_ ="); 

    double eta = max(delta_  /den_);
    double a   = max(delta1_ /den_);
    double be  = max(delta2_ /den_);
    double c   = max(delta3_ /den_);
    double de  = max(delta12_ /den_);
    double e   = max(delta13_ /den_);
    double f   = max(delta23_ /den_);
    double d1  = max(Delta_.col(0) / den_ );
    double d2  = max(Delta_.col(1) / den_ );
    double d3  = max(Delta_.col(2) / den_ );

    mat Lambda_; Lambda_.zeros(3,3);
    Lambda_ << a << de << e << endr
            << 0 << be << f << endr
            << 0 << 0  << c << endr;
    vec k_ = {d1,d2,d3};

    cout << "eta = " << eta << endl;
    //cout << Lambda_ << endl;
    //cout << (vec){d1,d2} << endl;
    Lambda_.print("Lambda_ =");
    k_.print("k_ =");

    mat iI_Delta_ = inv(eye(3,3) - Delta_);
    vec eta_  = iI_Delta_*(delta_ + 20.0*ones(3,1) ); 
    vec eta1_  = iI_Delta_*delta1_; 
    vec eta2_  = iI_Delta_*delta2_; 
    vec eta3_  = iI_Delta_*delta3_; 
    vec eta12_ = iI_Delta_*delta12_; 
    vec eta13_ = iI_Delta_*delta13_; 
    vec eta23_ = iI_Delta_*delta23_; 
    mat Gamma_ = iI_Delta_*Delta_;
    iI_Delta_.print("iI_Delta_ =");
    eta_.print("eta_ =");
    eta1_.print("eta1_ =");
    eta2_.print("eta2_ =");
    eta3_.print("eta3_ =");
    eta12_.print("eta12_ =");
    eta13_.print("eta13_ =");
    eta23_.print("eta23_ =");
    Gamma_.print("Gamma_ =");

    cube Lambda__; Lambda__.zeros(3,3,3);
    Lambda__.slice(0) << eta1_(0) << eta12_(0) << eta13_(0) << endr
                      << 0        << eta2_(0)  << eta23_(0) << endr
                      << 0        << 0         << eta3_(0)  << endr;
    Lambda__.slice(1) << eta1_(1) << eta12_(1) << eta13_(1) << endr
                      << 0        << eta2_(1)  << eta23_(1) << endr
                      << 0        << 0         << eta3_(1)  << endr;
    Lambda__.slice(2) << eta1_(2) << eta12_(2) << eta13_(2) << endr
                      << 0        << eta2_(2)  << eta23_(2) << endr
                      << 0        << 0         << eta3_(2)  << endr;

    Lambda__.print("Lambda__ =");

 
    return 0;
}