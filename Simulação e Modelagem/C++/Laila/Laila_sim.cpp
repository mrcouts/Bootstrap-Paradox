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
    double r = 0.37;
    double x0 = 0.0;
    double y0 = 0.0;
    double w = 1/r;
    return {x0+r*cos(w*t), y0+r*sin(w*t),0.480};
    //return {0,0,0.480};
}

vec dr_(double t){
    double r = 0.37;
    double w = 1/r;
    return {-w*r*sin(w*t), w*r*cos(w*t),0};
    //return {0,0,0};
}

vec d2r_(double t){
    double r = 0.37;
    double w = 1/r;
    return {-w*w*r*cos(w*t), -w*w*r*sin(w*t),0};
    //return {0,0,0};
}

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

    arma_rng::set_seed_random();
    vec sigma_ = {-0.15,0.15,0.15,0.15,0.15,0.15,-0.15,0.15,0.15,0.15,0.15,0.15,0.15};
    vec coef_  = ones(13) + sigma_;

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

    vec _l_ = {0.305, 0.045, 0.045 , 0.372, 0.0125, 0.000};
    vec _lg_= {coef_(0)*lg1  , 0.025, 0.0116, coef_(1)*lg4  , 0.0042, 0.000};
    vec _m_ = {coef_(2)*m1   , 0.033, 0.034 , coef_(3)*m4   , 0.023 , 0.000};

    cube _I__; _I__.zeros(3,3,6);
    _I__.slice(0) << coef_(4)*0.99*J1 << 0             << 0       << endr
                  << 0                << coef_(4)*J1   << 0       << endr
                  << 0                << 0             << 12.7e-6 << endr;

    _I__.slice(1) << 8.0e-6 << 0      << 0      << endr
                  << 0      << 1.5e-6 << 0      << endr
                  << 0      << 0      << 8.2e-6 << endr;

    _I__.slice(2) << 7.9e-6 << 0      << 0      << endr
                  << 0      << 8.2e-6 << 0      << endr
                  << 0      << 0      << 1.2e-6 << endr;

    _I__.slice(3) << coef_(5)*J4 << 0      << 0           << endr
                  << 0           << 1.7e-6 << 0           << endr
                  << 0           << 0      << coef_(5)*J4 << endr;

    _I__.slice(4) << 1.6e-6 << 0      << 0      << endr
                  << 0      << 1.3e-6 << 0      << endr
                  << 0      << 0      << 1.2e-6 << endr;

    _I__.slice(5) << 0 << 0 << 0 << endr
                  << 0 << 0 << 0 << endr
                  << 0 << 0 << 0 << endr;

    vec _l2_ = {0.000, 0.380, 0.000, 0.200};
    vec _lg2_= {0.000, coef_(6)*lg22 , 0.000, coef_(7)*lg42 };
    vec _m2_ = {coef_(8)*m12  , coef_(9)*m22  , 0.000, coef_(10)*m42  };

    cube _I2__; _I2__.zeros(3,3,4);
    _I2__.slice(0) << 0  << 0  << 0 << endr
                   << 0  << 0  << 0 << endr
                   << 0  << 0  << 0 << endr;

    _I2__.slice(1) << 0 << 0 << 0             << endr
                   << 0 << 0 << 0             << endr
                   << 0 << 0 << coef_(11)*J22 << endr;

    _I2__.slice(2) << 0  << 0  << 0 << endr
                   << 0  << 0  << 0 << endr
                   << 0  << 0  << 0 << endr;

    _I2__.slice(3) << 0  << 0  << 0 << endr
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

    Mecanismo _P = Mecanismo(3); _P.dy->Mh_ = coef_(12)*m0*eye(3,3);
    Serial __6R1 = Serial(6, _l_, _lg_ , _m_, _I__ , {0, 9.8,0}, &fDH_6R);
    Serial __6R2 = Serial(6, _l_, _lg_ , _m_, _I__ , {0, 9.8,0}, &fDH_6R);
    Serial _PRRP = Serial(4, _l2_, _lg2_ , _m2_, _I2__ , {9.8,0, 0}, &fDH_PRRP);
    Serial **_R_ = new Serial* [3];
    _R_[0] = &__6R1;
    _R_[1] = &__6R2;
    _R_[2] = &_PRRP;

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
    Parallel _Robot = Parallel(3, &_P, _R_, 3, {3,9,15}, D_, E_, F_, f_, P_, Q_, S_);
    Reference RefObj = Reference(1.0, {0.0, 0.0, 0.480}, {0.37, 0.0, 0.480});

    double lambda = 50.0;

    vec eta_ = {32.7299, 50.8718, 63.2811};
    cube Lambda__; Lambda__.zeros(3,3,3);
    Lambda__.slice(0) << 4.0933e+02 << 2.3252e+02 << 3.8963e+02 << endr
                      << 0          << 3.3300e+01 << 1.1099e+02 << endr
                      << 0          << 0          << 9.3153e+01 << endr;
    Lambda__.slice(1) << 9.5456e+02 << 5.4236e+02 << 9.0941e+02 << endr
                      << 0          << 7.7657e+01 << 2.5886e+02 << endr
                      << 0          << 0          << 2.1734e+02 << endr;
    Lambda__.slice(2) << 1.4451e+03 << 8.2065e+02 << 1.3759e+03 << endr
                      << 0          << 1.1743e+02 << 3.9178e+02 << endr
                      << 0          << 0          << 3.2897e+02 << endr;
    mat Gamma_; Gamma_.zeros(3,3);
    Gamma_ << 0.3018 << 0.1876 << 0.1021 << endr
           << 0.4071 << 0.8102 << 0.1664 << endr
           << 0.8604 << 0.7925 << 0.3145 << endr;
    //SMCLaw SMC = SMCLaw(3, lambda, 10.0, zeros(3,3), zeros(3), 100.0, &RefObj, &Robot);
    //SMCLaw SMC = SMCLaw(3, lambda, 10.0, 20.0, &r_, &dr_, &d2r_, &Robot);
    SMCLaw SMC = SMCLaw(3, lambda, eta_, Lambda__, Gamma_, 20.0, &r_, &dr_, &d2r_, &Robot);
    //SMCLaw SMC = SMCLaw(3, lambda, eta_, Lambda__, Gamma_, 20.0, &RefObj, &Robot);
    //FLControlLaw FL = FLControlLaw(3, lambda*lambda, 2*lambda, &RefObj, &Robot);
    //FLControlLaw FL = FLControlLaw(3, 400.0, 40.0, &r_, &dr_, &d2r_, &Robot);
    Acceleration AC = Acceleration(19, &_Robot, &SMC);
    //Acceleration AC = Acceleration(19, &Robot, &FL);


    //vec q0_ = {0.0, 0.0, 0.480,
    //           9.6321e-01, 0, 1.7529, 0, -1.1452, 0,
    //           9.6321e-01, 0, 1.7529, 0, -1.1452, 0,
    //           0, PI/2, 0, 9.6100e-02};
    vec q0_ = {3.7000e-01,
               1.0020e-17,
               4.8000e-01,
               1.6348e+00,
               1.1851e+00,
               1.2787e+00,
              -6.1400e-02,
              -1.3032e+00,
              -1.1681e+00,
               1.6348e+00,
              -1.1851e+00,
               1.2787e+00,
               6.1400e-02,
              -1.3032e+00,
               1.1681e+00,
               3.7000e-01,
               1.5708e+00,
              -1.2748e-16,
               9.6100e-02};
    vec x0_ = join_vert(q0_, zeros(19,1));

    /*
    Acceleration AC = Acceleration(16, &Robot, &RefObj);
    vec q0_ = {9.6321e-01, 0, 1.7529, 0, -1.1452, 0,
               9.6321e-01, 0, 1.7529, 0, -1.1452, 0,
               0, PI/2, 0, 9.6100e-02};
    vec x0_ = join_vert(q0_, zeros(16,1));
    */
    
    RK rk = RK("RK8", &AC);
    rk.Doit(0.000025, 4.65, x0_);
    double t;
    vec ref_; ref_.zeros(3);
    vec dref_; dref_.zeros(3);
    for(uint i = 0; i< rk.t_.n_rows; i++){
        t = rk.t_(i);
        cout << t << "; " << rk.u__(0,0,i)  << "; " << rk.u__(1,0,i) << "; " << rk.u__(2,0,i) << "; ";
        ref_ = r_(t);
        dref_ = dr_(t);
        //RefObj.Doit(rk.t_(i));
        //ref_ = RefObj.r_;
        //dref_ = RefObj.dr_;
        cout << ref_(0) - rk.y__(0,0,i)  << "; " << ref_(1) - rk.y__(1,0,i) << "; " << ref_(2) - rk.y__(2,0,i) << "; ";
        cout << -( dref_(0) - rk.y__(19,0,i) + lambda*(ref_(0) - rk.y__(0,0,i)) ) << "; " << -( dref_(1) - rk.y__(20,0,i) + lambda*(ref_(1) - rk.y__(1,0,i))) << "; " << -( dref_(2) - rk.y__(21,0,i) + lambda*(ref_(2) - rk.y__(2,0,i))) << "; " << endl;
    }
    //for(uint i = 0; i< rk.t_.n_rows; i++) cout << rk.t_(i) << "; " << rk.u__(0,0,i)  << "; " << rk.u__(1,0,i) << "; " << rk.u__(2,0,i) << "; " << endl;
    return 0;
}