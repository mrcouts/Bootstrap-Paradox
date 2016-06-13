#include "Acceleration.h"

Acceleration::Acceleration(int dof, Serial *R, FLControlLaw *u){
	this->dof = dof;
	this->R = R;
	this->u = u; }

Acceleration::~Acceleration(){}

field<vec> Acceleration::Doit(double t, vec q0_, vec q1_){
	Dy *dy = R->Doit(q0_, q1_);
	field<vec> F(2);
	F(1) = u->Doit(t, q0_, q1_);
	F(0) = solve(dy->Mh_, F(1) - dy->vh_ - dy->gh_ );
	return F; }

vec Acceleration::f_(double t, vec y_){
	return join_vert(y_(span(dof,2*dof-1)), Doit(t, y_(span(0,dof-1)), y_(span(dof,2*dof-1)))(0) ) ; }

field<vec> Acceleration::f2_(double t, vec y_){
	field<vec> F(2);
	F = Doit(t, y_(span(0,dof-1)), y_(span(dof,2*dof-1)));
	field<vec> F2(2);
	F2(0) = join_vert(y_(span(dof,2*dof-1)), F(0) );
	F2(1) = F(1); 
	return F2; }