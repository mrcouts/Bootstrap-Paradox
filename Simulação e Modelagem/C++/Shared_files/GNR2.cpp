#include "GNR2.h"
#include <iostream>

GNR2::GNR2(string method, vec x0_, vec (*f_)(vec), mat (*J_)(vec), double tol, uint nmax):GNR(method, x0_, f_, J_, tol, nmax){
	rk = new RK(method, this);
}

GNR2::~GNR2(){
	delete rk;
}

void GNR2::Doit(){
	for(uint i=0; i<nmax; i++){
		if(convergiu) return;
		rk->Doit(1,1, x0_);
		x_ = rk->y__.slice(1);
		//x_ = x0_ - solve(J_(x0_), res_);
		res_ = f_(x_);
		if( abs(norm(res_, "inf")) <  abs(tol) ){
			convergiu = true;
			n = i+1;
		}
		else
			x0_ = x_;
	}
}