#include <iostream>
#include <iomanip>
#include <string>
#include "RK.h"

vec f_(double t, vec y_){
	return {y_(1), -20*y_(1) - 100*y_(0)}; }

int main(void){
    RK RK6 = RK("RK6");
    RK6.Doit(0.1, 1, {100.0, 0}, &f_);
    cout << RK6.y__ << endl;
    return 0; }