#include "Acceleration.h"

Acceleration::Acceleration(int dof, Serial *R, FLControlLaw *u){
	this->dof = dof;
	this->R = R;
	this->u = u; }

Acceleration::~Acceleration(){}

vec Acceleration::Doit(double t, vec q0_, vec q1_){
	Dy *dy = R->Doit(q0_, q1_);
	return solve(dy->Mh_, u->Doit(t, q0_, q1_) - dy->vh_ - dy->gh_ ); }

vec Acceleration::f_(double t, vec q0_, vec q1_){
	return join_vert(q1_, Doit(t, q0_, q1_) ) ; }