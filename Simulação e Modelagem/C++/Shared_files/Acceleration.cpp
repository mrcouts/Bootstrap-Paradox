#include "Acceleration.h"

Acceleration::Acceleration(int dof, Serial *R, FLControlLaw *u){
	this->dof = dof;
	this->R = R;
	this->u = u;
	caso = 1; }

Acceleration::Acceleration(int dof, Parallel *R2, FLControlLaw *u){
	this->dof = dof;
	this->R2 = R2;
	this->u = u;
	caso = 2; }

Acceleration::~Acceleration(){}

field<vec> Acceleration::Doit(double t, vec q0_, vec q1_){
	Dy *dy;
	field<vec> F(2);
	switch(caso){
		case 1: 
		    dy = R->Doit(q0_, q1_);
		    F(1) = u->Doit(t, q0_, q1_);
		    F(0) = solve(dy->Mh_, F(1) - dy->vh_ - dy->gh_ );
		    break;
		case 2:
		    dy = R2->Doit(q0_, q1_);
		    F(1) = u->Doit(t, q0_ , q1_);
		    double lambda = 100.0;
		    if(R2->caso == 1) F(0) = solve(join_vert( R2->C_.t()*R2->M_, R2->A_ ), join_vert( R2->Z_.t()*F(1) -R2->C_.t()*(R2->v_ + R2->g_), R2->b_ - 2*lambda*R2->A_*q1_ - lambda*lambda*R2->_q_ ));
		    if(R2->caso == 2) F(0) = solve(join_vert( R2->C_.t()*R2->M_, R2->A_ ), join_vert( R2->Z_.t()*F(1) -R2->C_.t()*(R2->v_ + R2->g_), R2->b_ - lambda*join_diag(2*eye(R2->_q_.n_rows,R2->_q_.n_rows), eye(R2->b_.n_rows - R2->_q_.n_rows,R2->b_.n_rows - R2->_q_.n_rows))*R2->A_*q1_ - lambda*lambda*join_vert(R2->_q_, zeros(R2->b_.n_rows - R2->_q_.n_rows) ) ));
		    break; }
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