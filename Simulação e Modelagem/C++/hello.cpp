#include <iostream>
#include <iomanip>
#include "SomeUtilities.h"
#include "Serial.h"

#define PI 3.14159265358979323846264338327950288

mat fDH(vec q0_, vec l_, vec lg_){
    mat M;
    M << l_(0) << 0 << 0 << q0_(0) << -l_(0) + lg_(0) << 0 << 0 << true << endr
      << l_(1) << 0 << 0 << q0_(1) << -l_(1) + lg_(1) << 0 << 0 << true << endr;
    return M;}

//Simple C++ program
int main(void)
{
	cout << setprecision(16);
    cout << "Hello! This is a C++ program." << endl;
    cout << H_d(1, 0.1, 2, 0.2) << endl;

    cube I__; I__.zeros(3,3,2);
    I__.slice(0) << 0 << 0      << 0      << endr
                 << 0 << 0.0001 << 0      << endr
                 << 0 << 0      << 0.0001 << endr;

    I__.slice(1) << 0 << 0      << 0      << endr
                 << 0 << 0.0001 << 0      << endr
                 << 0 << 0      << 0.0001 << endr;

    //bool array[] = {true, true};
    Serial RR = Serial(2, {0.1, 0.1}, {0.05, 0.05},{0.1, 0.1}, I__ , {0, -9.8,0}, &fDH);

    //double t = 0;
    //vec q0_;
    //vec q1_;
    //vec q2_;
    //double w1 = 10;
    //double w2 = 15;
    //for(int i = 0; i<10000; i++){
    //	q0_ = {sin(w1*t),sin(w2*t)};
    //	q1_ = {w1*cos(w1*t),w2*cos(w2*t)};
    //	q2_ = {-w1*w1*sin(w1*t), -w2*w2*sin(w2*t)};
    //	RR.Doit(q0_,q1_);
    //	cout << RR.Mh_*q2_ + RR.vh_ + RR.gh_ << endl;
    //	//cout << RR.Mh_ << endl;
    //	//cout << RR.vh_ << endl;
    //	//cout << RR.gh_ << endl;
    //	t += 0.01;
    //}

    RR.Doit({PI/6,PI/4},{20.0,30.0});
    cout << RR.Mh_ << endl;
    cout << RR.vh_ << endl;
    cout << RR.gh_ << endl;
    RR.Mh_.raw_print(cout, "Mh_=");
    RR.vh_.raw_print(cout, "vh_=");
    RR.gh_.raw_print(cout, "gh_=");


    I__.clear();

    //cout << RR.l_  << endl;
    //cout << RR.lg_ << endl;
    //cout << RR.Hr__ << endl;
    //cout << RR.H__ << endl;
    //cout << RR.z__ << endl;
    //cout << RR.o__ << endl;
    //cout << RR.og__ << endl;
    //cout << RR.fDH(RR.l_, RR.l_, RR.lg_) << endl;
    //cout << RR.Jv__ << endl;
    //cout << RR.Jw__ << endl;
    //cout << RR.Jv_n_ << endl;
    //cout << RR.Jw_n_ << endl;
    //cout << RR.m_ << endl;
    //cout << RR.I__ << endl;
    //cout << RR.Jw2__ << endl;
    //cout << RR.M_ << endl;
    //cout << RR.w_rel__ << endl;
    //cout << RR.w_arr__ << endl;
    //vec v_rel = RR.Jv__(span(0,2),span(1,1),span(1,1))*30.0;
    //cout << v_rel << endl;
    //cout << 2*cross(RR.w_arr__.slice(1), v_rel);
    //cout << RR.dw_co__ << endl;
    //cout << RR.a_cen__ << endl;
    //cout << RR.a_cor__ << endl;
    //cout << RR.a_co__ << endl;
    //cout << RR.dw_co2__ << endl;


    return 0;
}