#include "GNR.h"
#include <iostream>

GNR::GNR(vec (*f_)(vec), mat (*J_)(vec), double tol, uint nmax){
	this->f_ = f_;
	this->J_ = J_;
	this->tol = tol;
	this->nmax = nmax;


	convergiu = false;
	n = nmax;
	caso = 1;
}

GNR::GNR(Parallel *R, double tol, uint nmax){
	this->R = R;
	this->tol = tol;
	this->nmax = nmax;

	convergiu = false;
	n = nmax;
	caso = 2;
}

GNR::~GNR(){}

vec GNR::g_(vec y_){
	switch(caso){
		case 1: return -solve(J_(y_), f_(x0_) );
		case 2: R->DoitKin( join_vert(R->P->q0_, y_), zeros(R->nq) );
		        return -solve(R->Ao1_ , R->_q_ );
	}
}

void GNR::Doit(vec x0_){
	this->x0_ = x0_;
	x_ = x0_;
	res_ = f_(x0_);
	if( abs(norm(res_, "inf" )) <  abs(tol) )
		convergiu = true;

	for(uint i=0; i<nmax; i++){
		if(convergiu) return;
		x_ = this->x0_ - solve(J_(this->x0_), res_);
		res_ = f_(x_);
		if( abs(norm(res_, "inf")) <  abs(tol) ){
			convergiu = true;
			n = i+1;
		}
		else
			this->x0_ = x_;
	}
}
