#include <iostream>
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

    double theta = PI/4;
    for(int i = 0; i<10000; i++){
    	RR.Doit({theta, theta});
    	cout << RR.M_ << endl;
    	cout << RR.gh_ << endl;
    	theta += 0.03;
    }

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

    return 0;
}