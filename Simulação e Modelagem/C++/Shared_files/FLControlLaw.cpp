#include "FLControlLaw.h"

ControlLaw::ControlLaw(uint dof, double lambda, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Serial *R){
	this->dof = dof;
	this->lambda = lambda;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->R = R;
	caso = 1;
	smc = false;
	RefObjFlag = false;
}

ControlLaw::ControlLaw(uint dof, double lambda, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Dy* (*dy_comp)(vec, vec)){
	this->dof = dof;
	this->lambda = lambda;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->dy_comp = dy_comp;
	caso = 2;
	smc = false;
	RefObjFlag = false;
}

ControlLaw::ControlLaw(uint dof, double lambda, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2){
	this->dof = dof;
	this->lambda = lambda;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->R2 = R2;
	caso = 3;
	smc = false;
	RefObjFlag = false;
}

ControlLaw::ControlLaw(uint dof, double lambda, Reference *RefObj, Serial *R){
	this->dof = dof;
	this->lambda = lambda;
	this->RefObj = RefObj;
	this->R = R;
	caso = 1;
	smc = false;
	RefObjFlag = true;
}

ControlLaw::ControlLaw(uint dof, double lambda, Reference *RefObj, Dy* (*dy_comp)(vec, vec)){
	this->dof = dof;
	this->lambda = lambda;
	this->RefObj = RefObj;
	this->dy_comp = dy_comp;
	caso = 2;
	smc = false;
	RefObjFlag = true;
}

ControlLaw::ControlLaw(uint dof, double lambda, Reference *RefObj, Parallel *R2){
	this->dof = dof;
	this->lambda = lambda;
	this->RefObj = RefObj;
	this->R2 = R2;
	caso = 3;
	smc = false;
	RefObjFlag = true;
}

ControlLaw::ControlLaw(uint dof, double lambda, double K, double phi, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Serial *R){
	this->dof = dof;
	this->lambda = lambda;
	this->K = K;
	this->phi = phi;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->R = R;
	caso = 1;
	smc = true;
	RefObjFlag = false;
}

ControlLaw::ControlLaw(uint dof, double lambda, double K, double phi, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Dy* (*dy_comp)(vec, vec)){
	this->dof = dof;
	this->lambda = lambda;
	this->K = K;
	this->phi = phi;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->dy_comp = dy_comp;
	caso = 2;
	smc = true;
	RefObjFlag = false;
}

ControlLaw::ControlLaw(uint dof, double lambda, double K, double phi, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2){
	this->dof = dof;
	this->lambda = lambda;
	this->K = K;
	this->phi = phi;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->R2 = R2;
	caso = 3;
	smc = true;
	RefObjFlag = false;
}

ControlLaw::ControlLaw(uint dof, double lambda, double K, double phi, Reference *RefObj, Serial *R){
	this->dof = dof;
	this->lambda = lambda;
	this->K = K;
	this->phi = phi;
	this->RefObj = RefObj;
	this->R = R;
	caso = 1;
	smc = true;
	RefObjFlag = true;
}

ControlLaw::ControlLaw(uint dof, double lambda, double K, double phi, Reference *RefObj, Dy* (*dy_comp)(vec, vec)){
	this->dof = dof;
	this->lambda = lambda;
	this->K = K;
	this->phi = phi;
	this->RefObj = RefObj;
	this->dy_comp = dy_comp;
	caso = 2;
	smc = true;
	RefObjFlag = true;
}

ControlLaw::ControlLaw(uint dof, double lambda, double K, double phi, Reference *RefObj, Parallel *R2){
	this->dof = dof;
	this->lambda = lambda;
	this->K = K;
	this->phi = phi;
	this->RefObj = RefObj;
	this->R2 = R2;
	caso = 3;
	smc = true;
	RefObjFlag = true;
}

ControlLaw::~ControlLaw(){}

vec ControlLaw::Doit(double t, vec q0_, vec q1_){
	Dy *dy;

	if(smc)	{
		if(RefObjFlag){
			RefObj->Doit(t);
				switch (caso){
				case 1: dy = R->Doit(q0_, q1_); break;
				case 2: dy = dy_comp(q0_, q1_); break;
				case 3: dy = R2->Doit(q0_, q1_);
				        return dy->vh_ + dy->gh_ + dy->Mh_*(RefObj->d2r_ + lambda*(RefObj->dr_ - R2->P->q1_ ) + K*tanh(( RefObj->dr_ - R2->P->q1_  + lambda*(RefObj->r_ - R2->P->q0_)  )/phi) );
			}
			return             dy->vh_ + dy->gh_ + dy->Mh_*(RefObj->d2r_ + lambda*(RefObj->dr_ - q1_)         + K*tanh(( RefObj->dr_ - q1_         + lambda*(RefObj->r_ - q0_       )  )/phi) );
		}
		else{
			switch (caso){
				case 1: dy = R->Doit(q0_, q1_); break;
				case 2: dy = dy_comp(q0_, q1_); break;
				case 3: dy = R2->Doit(q0_, q1_);
				        return dy->vh_ + dy->gh_ + dy->Mh_*(d2r_(t) + lambda*(dr_(t) - R2->P->q1_ ) + K*tanh(( dr_(t) - R2->P->q1_ + lambda*(r_(t) - R2->P->q0_)  )/phi) );
			}
			return             dy->vh_ + dy->gh_ + dy->Mh_*(d2r_(t) + lambda*(dr_(t) - q1_)         + K*tanh(( dr_(t) - q1_        + lambda*(r_(t) - q0_       )  )/phi) );
		}
	}
	else{
		if(RefObjFlag){
			RefObj->Doit(t);
				switch (caso){
				case 1: dy = R->Doit(q0_, q1_); break;
				case 2: dy = dy_comp(q0_, q1_); break;
				case 3: dy = R2->Doit(q0_, q1_);
				        return dy->vh_ + dy->gh_ + dy->Mh_*(RefObj->d2r_ + 2*lambda*(RefObj->dr_ - R2->P->q1_ ) + lambda*lambda*(RefObj->r_ - R2->P->q0_  ) );
			}
			return             dy->vh_ + dy->gh_ + dy->Mh_*(RefObj->d2r_ + 2*lambda*(RefObj->dr_ - q1_)         + lambda*lambda*(RefObj->r_ - q0_) );
		}
		else{
			switch (caso){
				case 1: dy = R->Doit(q0_, q1_); break;
				case 2: dy = dy_comp(q0_, q1_); break;
				case 3: dy = R2->Doit(q0_, q1_);
				        return dy->vh_ + dy->gh_ + dy->Mh_*(d2r_(t) + 2*lambda*(dr_(t) - R2->P->q1_ ) + lambda*lambda*(r_(t) - R2->P->q0_  ) );
			}
			return dy->vh_ + dy->gh_ + dy->Mh_*(d2r_(t) + 2*lambda*(dr_(t) - q1_) + lambda*lambda*(r_(t) - q0_) );
		}
	}
}

FLControlLaw::FLControlLaw(uint dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Serial *R){
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

FLControlLaw::FLControlLaw(uint dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Dy* (*dy_comp)(vec, vec)){
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

FLControlLaw::FLControlLaw(uint dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2){
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

FLControlLaw::FLControlLaw(uint dof, double Kp, double Kv, Reference *RefObj, Serial *R){
	this->dof = dof;
	this->Kp = Kp;
	this->Kv = Kv;
	this->RefObj = RefObj;
	this->R = R;
	caso = 1;
	RefObjFlag = true;
}

FLControlLaw::FLControlLaw(uint dof, double Kp, double Kv, Reference *RefObj, Dy* (*dy_comp)(vec, vec)){
	this->dof = dof;
	this->Kp = Kp;
	this->Kv = Kv;
	this->RefObj = RefObj;
	this->dy_comp = dy_comp;
	caso = 2;
	RefObjFlag = true;
}

FLControlLaw::FLControlLaw(uint dof, double Kp, double Kv, Reference *RefObj, Parallel *R2){
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

SMCLaw::SMCLaw(uint dof, double Kp, double eta, double n, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2){
	this->dof = dof;
	this->Kp = Kp;
	this->eta = eta;
	this->n = n;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->R2 = R2;

	s_.zeros(dof);
	sigma_.zeros(dof);
	k.zeros(1,1);
	k_.zeros(dof);

	RefObjFlag = false;
	caso = 1;
}

SMCLaw::SMCLaw(uint dof, double Kp, double eta, double n, Reference *RefObj, Parallel *R2){
	this->dof = dof;
	this->Kp = Kp;
	this->eta = eta;
	this->n = n;
	this->RefObj = RefObj;
	this->R2 = R2;

	s_.zeros(dof);
	sigma_.zeros(dof);
	k.zeros(1,1);
	k_.zeros(dof);

	RefObjFlag = true;
	caso = 1;
}

SMCLaw::SMCLaw(uint dof, double Kp, double eta, mat Lambda_, vec gamma_, double n, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2){
	this->dof = dof;
	this->Kp = Kp;
	this->eta = eta;
	this->Lambda_ = Lambda_;
	this->gamma_ = gamma_;
	this->n = n;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->R2 = R2;

	s_.zeros(dof);
	sigma_.zeros(dof);
	k.zeros(1,1);
	k_.zeros(dof);

	RefObjFlag = false;
	caso = 2;
}

SMCLaw::SMCLaw(uint dof, double Kp, double eta, mat Lambda_, vec gamma_, double n, Reference *RefObj, Parallel *R2){
	this->dof = dof;
	this->Kp = Kp;
	this->eta = eta;
	this->Lambda_ = Lambda_;
	this->gamma_ = gamma_;
	this->n = n;
	this->RefObj = RefObj;
	this->R2 = R2;

	s_.zeros(dof);
	sigma_.zeros(dof);
	k.zeros(1,1);
	k_.zeros(dof);

	RefObjFlag = true;
	caso = 2;
}

SMCLaw::SMCLaw(uint dof, double Kp, vec eta_, cube Lambda__, mat Gamma_, double n, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Parallel *R2){
	this->dof = dof;
	this->Kp = Kp;
	this->eta_ = eta_;
	this->Lambda__ = Lambda__;
	this->Gamma_ = Gamma_;
	this->n = n;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->R2 = R2;

	s_.zeros(dof);
	sigma_.zeros(dof);
	k.zeros(1,1);
	k_.zeros(dof);

	RefObjFlag = false;
	caso = 3;
}

SMCLaw::SMCLaw(uint dof, double Kp, vec eta_, cube Lambda__, mat Gamma_, double n, Reference *RefObj, Parallel *R2){
	this->dof = dof;
	this->Kp = Kp;
	this->eta_ = eta_;
	this->Lambda__ = Lambda__;
	this->Gamma_ = Gamma_;
	this->n = n;
	this->RefObj = RefObj;
	this->R2 = R2;

	s_.zeros(dof);
	sigma_.zeros(dof);
	k.zeros(1,1);
	k_.zeros(dof);

	RefObjFlag = true;
	caso = 3;
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
	switch(caso){
		case 1: k_ = eta*ones(dof);
		        break;
		case 2: k = eta + abs(R2->P->q1_).t()*Lambda_*abs(R2->P->q1_) + gamma_.t()*abs(sigma_);
		        k_ = k(0,0)*ones(dof);
		        break;
		case 3: k_ = eta_ + Gamma_*abs(sigma_);
				for(uint i=0; i<dof; i++){
					k = abs(R2->P->q1_).t()*Lambda__.slice(i)*abs(R2->P->q1_);
					k_(i) += k(0,0);
				}
				break;
	}
	return  dy->vh_ + dy->gh_ + dy->Mh_*(sigma_ - k_ % tanh(n*s_) );
}