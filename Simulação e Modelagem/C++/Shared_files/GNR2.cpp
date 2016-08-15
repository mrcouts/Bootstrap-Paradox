#include "GNR2.h"
#include <iostream>

GNR2::GNR2(string method, vec (*f_)(vec), mat (*J_)(vec), double tol, uint nmax):GNR(f_, J_, tol, nmax){
	rk = new RK(method, this);
}

GNR2::GNR2(string method, Parallel *R, double tol, uint nmax):GNR(R, tol, nmax){
	rk = new RK(method, this);
}


GNR2::~GNR2(){
	delete rk;
}

void GNR2::Doit(vec x0_){
	this->x0_ = x0_;
	x_ = x0_;
	res_ = f_(x0_);
	if( abs(norm(res_, "inf" )) <  abs(tol) )
		convergiu = true;

	for(uint i=0; i<nmax; i++){
		if(convergiu) return;
		rk->Doit(1,1, this->x0_);
		x_ = rk->y__.slice(1);
		res_ = f_(x_);
		if( abs(norm(res_, "inf")) <  abs(tol) ){
			convergiu = true;
			n = i+1;
		}
		else
			this->x0_ = x_;
	}
}

void GNR2::Doit(vec x0_, vec q0h_){
	this->x0_ = x0_;
	R->P->q0_ = q0h_;
	x_ = x0_;
	R->DoitKin( join_vert(q0h_, x0_), zeros(R->nq) );
	res_ = R->_q_;
	if( abs(norm(res_, "inf" )) <  abs(tol) )
		convergiu = true;

	for(uint i=0; i<nmax; i++){
		if(convergiu) return;
		rk->Doit(1,1, this->x0_);
		x_ = rk->y__.slice(1);
		res_ = R->_q_;
		if( abs(norm(res_, "inf")) <  abs(tol) ){
			convergiu = true;
			n = i+1;
		}
		else
			this->x0_ = x_;
	}
}