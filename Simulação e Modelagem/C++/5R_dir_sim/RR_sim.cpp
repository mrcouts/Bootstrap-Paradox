#include <iostream>
#include <iomanip>
#include <string>
#include "SomeUtilities.h"
#include "Serial.h"
#include "RR.h"
#include "FLControlLaw.h"
#include "Acceleration.h"
#include "RK.h"

vec r_(double t){
    double w1 = 10;
    double w2 = 15;
    return {sin(w1*t),sin(w2*t)}; }

vec dr_(double t){
    double w1 = 10;
    double w2 = 15;
    return {w1*cos(w1*t),w2*cos(w2*t)}; }

vec d2r_(double t){
    double w1 = 10;
    double w2 = 15;
    return {-w1*w1*sin(w1*t), -w2*w2*sin(w2*t)}; }

Dy* (dy_comp)(vec q0_, vec q1_){
    Dy *dy;
    dy = new Dy(2);
    return dy; }

class _5R_{
public:
    _5R_(){
        M_.zeros(6,6);
        v_.zeros(6);
        g_.zeros(6);
        _q_.zeros(4);
        Ah_.zeros(4,2);
        Ao_.zeros(4,4);
        C_.zeros(6,2);
        Mh_.zeros(2,2);
        gh_.zeros(2);

        cube I__; I__.zeros(3,3,2);
        I__.slice(0) << 0 << 0        << 0        << endr
                     << 0 << 171.6e-6 << 0        << endr
                     << 0 << 0        << 171.6e-6 << endr;

        I__.slice(1) << 0 << 0        << 0        << endr
                     << 0 << 320.6e-6 << 0        << endr
                     << 0 << 0        << 320.6e-6 << endr;

        RR1 = new Serial(2, {0.12, 0.15}, {0.5*0.12, 0.5*0.15},{0.143, 0.171}, I__ , {0, -9.8,0}, &fDH_RR);
        RR2 = new Serial(2, {0.12, 0.15}, {0.5*0.12, 0.5*0.15},{0.143, 0.171}, I__ , {0, -9.8,0}, &fDH_RR);
        P = new Mecanismo(2); 
    }
    ~_5R_(){}
    void Doit(vec q0_, vec q1_){
        P->Doit(q0_(span(0,1)), q1_(span(0,1)));
        RR1->Doit(q0_(span(2,3)), q1_(span(2,3)));
        RR2->Doit(q0_(span(4,5)), q1_(span(4,5)));
        M_ = join_diag(0*P->dy->Mh_, join_diag(RR1->Mh_, RR2->Mh_) );
        cout << M_ << endl;
        v_ = join_vert(P->dy->vh_, join_vert(RR1->vh_, RR2->vh_) );
        g_ = join_vert(P->dy->gh_, join_vert(RR1->gh_, RR2->gh_) );
        mat H1 = Hy(0,0.05,0,0);
        mat R1 = Roty(0);
        vec w1 = H1*join_vert(RR1->o__.slice(2),ones(1));
        mat H2 = Hy(PI,-0.05,0,0);
        mat R2 = Roty(PI);
        vec w2 = H2*join_vert(RR2->o__.slice(2),ones(1));
        _q_ = join_vert(q0_(span(0,1)) - w1(span(0,1)), q0_(span(0,1)) - w2(span(0,1))) ;
        cout << _q_ << endl;
        Ah_ = join_vert(eye(2,2),eye(2,2));
        cout << Ah_ << endl;
        mat M1 = R1*RR1->Jv_n_;
        mat M2 = R2*RR2->Jv_n_;
        Ao_ = join_diag(-M1(span(0,1),span(0,1)), -M2(span(0,1),span(0,1)));
        cout << Ao_ << endl;
        C_ = join_vert( eye(2,2), -solve(Ao_,Ah_) );
        cout << C_ << endl;
        Mh_ = C_.t()*M_*C_;
        cout << Mh_ << endl;
        gh_ = C_.t()*g_;
        cout << gh_ << endl;    
    }
    mat M_;
    vec v_;
    vec g_;
    vec _q_;
    mat Ah_;
    mat Ao_;
    mat C_;
    mat Mh_;
    vec gh_;
    Serial *RR1;
    Serial *RR2;
    Mecanismo *P;
};

int main(void){
    int dof = 2;
    cube I__; I__.zeros(3,3,dof);
    I__.slice(0) << 0 << 0        << 0        << endr
                 << 0 << 171.6e-6 << 0        << endr
                 << 0 << 0        << 171.6e-6 << endr;

    I__.slice(1) << 0 << 0        << 0        << endr
                 << 0 << 320.6e-6 << 0        << endr
                 << 0 << 0        << 320.6e-6 << endr;

    Serial RR  = Serial(2, {0.12, 0.15}, {0.5*0.12, 0.5*0.15},{0.143, 0.171}, I__ , {0, -9.8,0}, &fDH_RR);

    FLControlLaw FL = FLControlLaw(dof, 100.0, 20.0, &r_, &dr_, &d2r_, &RR);
    Acceleration AC = Acceleration(dof, &RR, &FL);
    //vec u; u.zeros(dof);
    //u = FL.Doit(0, r_(0.1), dr_(0.1));
    //cout << u << endl;
    //vec v; v.zeros(2*dof);
    //v = AC.f_(0, join_vert(r_(0.1), dr_(0.1)) );
    //cout << v << endl;

    RK rk = RK("RK8", &AC);
    rk.Doit(0.001, 0.01, join_vert(r_(0.0), dr_(0.0)) );
    for(uint i = 0; i< rk.t_.n_rows; i++)
        cout << rk.t_(i) << "; " << rk.u__(0,0,i)  << "; " << rk.u__(1,0,i) << "; " << endl;
        //cout << rk.t_(i) << "; " << r_(rk.t_(i))(1) - rk.y__(1,0,i)  << "; " << dr_(rk.t_(i))(1) - rk.y__(3,0,i) << "; " << endl;

    //cout << RR.o__.slice(0) << endl;

    _5R_ Robot = _5R_();
    Robot.Doit({0.0, 0.20,PI/4+0.001,PI/4-0.001,PI/4-0.001,PI/4+0.001},{0,0,0,0,0,0});

    return 0; }