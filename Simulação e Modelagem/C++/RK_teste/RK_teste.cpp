#include <iostream>
#include <iomanip>
#include <string>
#include "RK.h"
#include "GNR2.h"

//vec f_(double t, vec y_){
//	//return {y_(1), -20*y_(1) - 100*y_(0)};
//	return {-(2-1)*(2-1)/(2*(y_(0)-1)) };
//}

vec f_(vec x_){
	return x_ % x_; 
}

mat J_(vec x_){
	return 2*diagmat(x_);
}

int main(void){
	GNR gnr = GNR("Lalala", {2}, &f_, &J_, 1e-7, 100);
    RK RK6 = RK("RK6", &gnr);
    //RK6.Doit(0.01, 1.5, {100.0, 0});
    RK6.Doit(1, 1, {2});
    for(uint i = 0; i< RK6.t_.n_rows; i++)
    	cout << RK6.t_(i) << "; " << RK6.y__(0,0,i)  << "; " << RK6.y__(0,0,i) << "; " << endl;
    	//cout << RK6.t_(i) << "; " << RK6.y__(0,0,i)  << "; " << RK6.y__(1,0,i) << "; " << endl;

    GNR2 gnr2 = GNR2("RK6", {2}, &f_, &J_, 1e-10, 100);
    gnr2.Doit();
    cout << gnr2.convergiu << endl;
    cout << gnr2.x_ << endl;
    cout << gnr2.res_ << endl;
    cout << gnr2.n << endl;
    return 0;
}