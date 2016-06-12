#include "FLControlLaw.h"

FLControlLaw::FLControlLaw(int dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Serial *R){
	this->dof = dof;
	this->Kp = Kp;
	this->Kv = Kv;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->R = R;
	caso = 1; }

FLControlLaw::FLControlLaw(int dof, double Kp, double Kv, vec (*r_)(double), vec (*dr_)(double), vec (*d2r_)(double), Dy* (*dy_comp)(vec, vec)){
	this->dof = dof;
	this->Kp = Kp;
	this->Kv = Kv;
	this->r_ = r_;
	this->dr_ = dr_;
	this->d2r_ = d2r_;
	this->dy_comp = dy_comp;
	caso = 2; }

FLControlLaw::~FLControlLaw(){}

vec FLControlLaw::Doit(double t, vec q0_, vec q1_){
	Dy *dy;
	switch (caso){
		case 1: dy = R->Doit(q0_, q1_); break;
		case 2: dy = dy_comp(q0_, q1_); break; }
	return dy->vh_ + dy->gh_ + dy->Mh_*(d2r_(t) + Kv*(dr_(t) - q1_) + Kp*(r_(t) - q0_) ); }