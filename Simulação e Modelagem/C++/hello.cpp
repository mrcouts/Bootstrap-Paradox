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

    //bool array[] = {true, true};
    Serial RR = Serial(2, {0.1, 0.1}, {0.05, 0.05}, &fDH);
    RR.Doit({PI/4, PI/4});
    //	RR.fDH = fDH;

    cout << RR.l_  << endl;
    cout << RR.lg_ << endl;
    cout << RR.Hr__ << endl;
    cout << RR.H__ << endl;
    cout << RR.z__ << endl;
    cout << RR.o__ << endl;
    cout << RR.og__ << endl;
    cout << RR.fDH(RR.l_, RR.l_, RR.lg_) << endl;

    cube C; C.zeros(3,3,3);
    cout << C << endl;
    vec v = {1,2,3};
    cout << v << endl;

    C( span(0,2), span(0,0), span(0,0) ) = v;
    cout << C << endl;

    cout << RR.Jv__ << endl;
    cout << RR.Jw__ << endl;


    return 0;
}