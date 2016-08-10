#include <iostream>
#include <iomanip>
#include <string>
#include "RK.h"

vec f_(double t, vec y_){
	//return {y_(1), -20*y_(1) - 100*y_(0)};
	return {-(2-1)*(2-1)/(2*(y_(0)-1)) };
}

int main(void){
    RK RK6 = RK("RK6", &f_);
    //RK6.Doit(0.01, 1.5, {100.0, 0});
    RK6.Doit(0.1, 1, {2});
    for(uint i = 0; i< RK6.t_.n_rows; i++)
    	cout << RK6.t_(i) << "; " << RK6.y__(0,0,i)  << "; " << RK6.y__(0,0,i) << "; " << endl;
    	//cout << RK6.t_(i) << "; " << RK6.y__(0,0,i)  << "; " << RK6.y__(1,0,i) << "; " << endl;
    return 0;
}