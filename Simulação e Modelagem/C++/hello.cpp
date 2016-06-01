#include <iostream>
#include "SomeUtilities.h"
#include "Serial.h"
#include <stdio.h>  /* for printf */

#define PI 3.14159265358979323846264338327950288

mat fDH(vec q0_, vec l_, vec lg_){
    mat M;
    M << l_(0) << 0 << 0 << q0_(0) << -l_(0) + lg_(0) << 0 << 0 << endr
      << l_(1) << 0 << 0 << q0_(1) << -l_(1) + lg_(1) << 0 << 0 << endr;
    return M;}

//Simple C++ program
int main(void)
{
    cout << "Hello! This is a C++ program." << endl;
    cout << H_d(1, 0.1, 2, 0.2) << endl;

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

    return 0;
}