#include "FLControlLaw.h"

FLControlLaw::FLControlLaw(int dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Serial *R){
	this->dof = dof;
	this->Kp = Kp;
	this->Kv = Kv;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->R = R;
	caso = 1;
	RefObjFlag = false;
}

FLControlLaw::FLControlLaw(int dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Dy* (*dy_comp)(vec, vec)){
	this->dof = dof;
	this->Kp = Kp;
	this->Kv = Kv;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->dy_comp = dy_comp;
	caso = 2;
	RefObjFlag = false;
}

FLControlLaw::FLControlLaw(int dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2){
	this->dof = dof;
	this->Kp = Kp;
	this->Kv = Kv;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->R2 = R2;
	caso = 3;
	RefObjFlag = false;
}

FLControlLaw::FLControlLaw(int dof, double Kp, double Kv, Reference *RefObj, Serial *R){
	this->dof = dof;
	this->Kp = Kp;
	this->Kv = Kv;
	this->RefObj = RefObj;
	this->R = R;
	caso = 1;
	RefObjFlag = true;
}

FLControlLaw::FLControlLaw(int dof, double Kp, double Kv, Reference *RefObj, Dy* (*dy_comp)(vec, vec)){
	this->dof = dof;
	this->Kp = Kp;
	this->Kv = Kv;
	this->RefObj = RefObj;
	this->dy_comp = dy_comp;
	caso = 2;
	RefObjFlag = true;
}

FLControlLaw::FLControlLaw(int dof, double Kp, double Kv, Reference *RefObj, Parallel *R2){
	this->dof = dof;
	this->Kp = Kp;
	this->Kv = Kv;
	this->RefObj = RefObj;
	this->R2 = R2;
	caso = 3;
	RefObjFlag = true;
}

FLControlLaw::~FLControlLaw(){}

vec FLControlLaw::Doit(double t, vec q0_, vec q1_){
	Dy *dy;
	if(RefObjFlag){
		RefObj->Doit(t);
			switch (caso){
			case 1: dy = R->Doit(q0_, q1_); break;
			case 2: dy = dy_comp(q0_, q1_); break;
			case 3: dy = R2->Doit(q0_, q1_);
			        return dy->vh_ + dy->gh_ + dy->Mh_*(RefObj->d2r_ + Kv*(RefObj->dr_ - R2->P->q1_ ) + Kp*(RefObj->r_ - R2->P->q0_  ) );
		}
		return dy->vh_ + dy->gh_ + dy->Mh_*(RefObj->d2r_ + Kv*(RefObj->dr_ - q1_) + Kp*(RefObj->r_ - q0_) );
	}
	else{
		switch (caso){
			case 1: dy = R->Doit(q0_, q1_); break;
			case 2: dy = dy_comp(q0_, q1_); break;
			case 3: dy = R2->Doit(q0_, q1_);
			        return dy->vh_ + dy->gh_ + dy->Mh_*(d2r_(t) + Kv*(dr_(t) - R2->P->q1_ ) + Kp*(r_(t) - R2->P->q0_  ) );
		}
		return dy->vh_ + dy->gh_ + dy->Mh_*(d2r_(t) + Kv*(dr_(t) - q1_) + Kp*(r_(t) - q0_) );
	}
}