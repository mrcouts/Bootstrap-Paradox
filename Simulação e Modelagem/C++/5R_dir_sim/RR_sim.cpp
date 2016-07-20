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
	int dof;
	int nq;
	mat M_;
    vec v_;
    vec g_;
    vec _q_;
    mat Ah_;
    mat Ao_;
    mat A_;
    vec b_;
    mat C_;
    mat Z_;
    mat Mh_;
    vec vh_;
    vec gh_;
    mat H_;
    mat h_;
    Mecanismo *P;
    Serial *RR1;
    Serial *RR2;
    _5R_(){
    	dof = 2;
    	nq = 6;
        M_.zeros(nq,nq);
        v_.zeros(nq);
        g_.zeros(nq);
        _q_.zeros(nq-dof);
        Ah_.zeros(nq-dof,dof);
        Ao_.zeros(nq-dof,nq-dof);
        A_.zeros(nq-dof,nq);
        b_.zeros(nq-dof);
        C_.zeros(nq,dof);
        Z_.zeros(dof,dof);
        Mh_.zeros(dof,dof);
        vh_.zeros(dof);
        gh_.zeros(dof);
        H_.zeros(dof,dof);
        h_.zeros(dof);

        cube I__; I__.zeros(3,3,2);
        I__.slice(0) << 0 << 0        << 0        << endr
                     << 0 << 171.6e-6 << 0        << endr
                     << 0 << 0        << 171.6e-6 << endr;

        I__.slice(1) << 0 << 0        << 0        << endr
                     << 0 << 320.6e-6 << 0        << endr
                     << 0 << 0        << 320.6e-6 << endr;

        RR1 = new Serial(2, {0.12, 0.15}, {0.5*0.12, 0.5*0.15},{0.143, 0.171}, I__ , {0, -9.8,0}, &fDH_RR);
        RR2 = new Serial(2, {0.12, 0.15}, {0.5*0.12, 0.5*0.15},{0.143, 0.171}, I__ , {0, -9.8,0}, &fDH_RR);
        P = new Mecanismo(dof); 
    }
    ~_5R_(){}
    void Doit(vec q0_, vec q1_){
    	vec q0h_ = q0_(span(0,1));
    	vec q1h_ = q1_(span(0,1));
    	vec q0h1_ = q0_(span(2,3));
    	vec q1h1_ = q1_(span(2,3));
    	vec q0h2_ = q0_(span(4,5));
    	vec q1h2_ = q1_(span(4,5));

        P->Doit(q0h_, q1h_);
        RR1->Doit(q0h1_, q1h1_);
        RR2->Doit(q0h2_, q1h2_);
        //M_, v_ e g_
        M_ = join_diag(P->dy->Mh_, join_diag(RR1->Mh_, RR2->Mh_) );
        v_ = join_vert(P->dy->vh_, join_vert(RR1->vh_, RR2->vh_) );
        g_ = join_vert(P->dy->gh_, join_vert(RR1->gh_, RR2->gh_) );
        // _q_, Ah_, Ao_ e A_
        mat H1 = Hy(0,0.05,0,0);
        mat R1 = Roty(0);
        vec w1 = H1*join_vert(RR1->o__.slice(RR1->dof),ones(1));
        mat H2 = Hy(PI,-0.05,0,0);
        mat R2 = Roty(PI);
        vec w2 = H2*join_vert(RR2->o__.slice(RR2->dof),ones(1));
        _q_ = join_vert(q0h_ - w1(span(0,1)), q0h_ - w2(span(0,1))) ;
        Ah_ = join_vert(eye(dof,dof),eye(dof,dof));
        mat M1 = R1*RR1->Jv_n_;
        mat M2 = R2*RR2->Jv_n_;
        Ao_ = join_diag(-M1(span(0,1),span(0,1)), -M2(span(0,1),span(0,1)));
        A_ = join_horiz(Ah_, Ao_);
        // b_
        vec v1 = R1*RR1->a_co_n_;
        vec v2 = R2*RR2->a_co_n_;
        b_ = join_vert( v1(span(0,1)), v2(span(0,1)) );
        // C_ e Z_
        C_ = join_vert( eye(dof,dof), -solve(Ao_,Ah_) );
        Z_ = join_vert(C_.row(2), C_.row(4));
        // Mh_, vh_ e gh_
        Mh_ = C_.t()*M_*C_;
        vh_ = C_.t()*( M_*join_vert(zeros(dof), solve(Ao_, b_) ) + v_ );
        gh_ = C_.t()*g_;
        // H_ e h_
        H_ = solve(Z_.t(), Mh_);
        h_ = solve(Z_.t(), vh_ + gh_);

        cout << Mh_ << endl;
        cout << vh_ << endl;
        cout << gh_ << endl;
        cout << C_ << endl;
        cout << Z_ << endl;
        cout << H_ << endl;
        cout << h_ << endl;
    }
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
    Robot.Doit({0.0, 0.20, 1.01377246756945, 1.41460649673445, 1.01377246756945, 1.41460649673445},{0,0,0,0,0,0});


    return 0;
}