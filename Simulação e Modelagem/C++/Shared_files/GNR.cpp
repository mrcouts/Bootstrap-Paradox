#include "GNR.h"
#include <iostream>

GNR::GNR(string method, vec x0_, vec (*f_)(vec), mat (*J_)(vec), double tol, uint nmax){
	this->method = method;
	this->x0_ = x0_;
	this->f_ = f_;
	this->J_ = J_;
	this->tol = tol;
	this->nmax = nmax;

	res_ = f_(x0_);
	if( abs(norm(res_, "inf" )) <  tol ) convergiu = true;
	else convergiu = false;
	x_ = x0_;
	n = nmax;
}

GNR::~GNR(){}

void GNR::Doit(){
	for(uint i=0; i<nmax; i++){
		if(convergiu) return;
		x_ = x0_ - solve(J_(x0_), res_);
		res_ = f_(x_);
		if( abs(norm(res_, "inf")) <  abs(tol) ){
			convergiu = true;
			n = i+1;
		}
		else
			x0_ = x_;
	}
}