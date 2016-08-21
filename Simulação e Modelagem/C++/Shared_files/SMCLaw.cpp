#include "SMCLaw.h"

SMCLaw::SMCLaw(int dof, double Kp, double eta, mat K_, vec k_, double n, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2){
	this->dof = dof;
	this->Kp = Kp;
	this->eta = eta;
	this->K_ = K_;
	this->k_ = k_;
	this->n = n;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->R2 = R2;
	RefObjFlag = false;

	s_.zeros(dof);
	sigma_.zeros(dof);
	k.zeros(1,1);
}

SMCLaw::SMCLaw(int dof, double Kp, double eta, mat K_, vec k_, double n, Reference *RefObj, Parallel *R2){
	this->dof = dof;
	this->Kp = Kp;
	this->eta = eta;
	this->K_ = K_;
	this->k_ = k_;
	this->n = n;
	this->RefObj = RefObj;
	this->R2 = R2;
	RefObjFlag = true;

	s_.zeros(dof);
	sigma_.zeros(dof);
	k.zeros(1,1);
}

SMCLaw::~SMCLaw(){}

vec SMCLaw::Doit(double t, vec q0_, vec q1_){
	Dy *dy;
	dy = R2->Doit(q0_, q1_);
	mat Aux_;
	if(RefObjFlag){
		RefObj->Doit(t);
		s_ = - ( (RefObj->dr_ - R2->P->q1_ ) + Kp*(RefObj->r_ - R2->P->q0_  ) );
		sigma_ = RefObj->d2r_ + Kp*(RefObj->dr_ - R2->P->q1_ );
	}
	else{
		s_ = - ( (dr_(t) - R2->P->q1_ ) + Kp*(r_(t) - R2->P->q0_  ) );
		sigma_ = d2r_(t) + Kp*(dr_(t) - R2->P->q1_ );
	}
	k = eta + abs(q1_).t()*K_*abs(q1_) + k_.t()*abs(sigma_);
	return  dy->vh_ + dy->gh_ + dy->Mh_*(sigma_ - k*tanh(n*s_) );
}