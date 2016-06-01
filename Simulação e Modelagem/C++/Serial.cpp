#include "Serial.h"

Serial::Serial(int dof, vec l_, vec lg_, mat (*fDH)(vec, vec, vec)){
	this->dof = dof;
	this->l_ = l_;
	this->lg_= lg_;
	this->Hr__.zeros(4,4,dof);
	this->H__.zeros(4,4,dof);
	this->z__.zeros(3,1,dof);
	this->o__.zeros(3,1,dof);
	this->og__.zeros(3,1,dof);
	this->fDH = fDH; }

Serial::~Serial(){
	l_.clear() ;
	lg_.clear();
	Hr__.clear();
	H__.clear();
	z__.clear();
	o__.clear();
	og__.clear(); }

void Serial::Doit(vec q0_){
	mat H;
	vec ogh_;
	for(int i = 0; i<dof; i++){
		H = this->fDH(q0_, l_, lg_);
		Hr__.slice(i) = H_d(H(i,0),H(i,1),H(i,2),H(i,3));
		H__.slice(i) = i==0 ? Hr__.slice(i): H__.slice(i-1)*Hr__.slice(i);
		z__.slice(i) = H__( span(0,2), span(2,2), span(i,i) );
		o__.slice(i) = H__( span(0,2), span(3,3), span(i,i) );
		ogh_ = H__.slice(i)*join_vert( H( span(i,i), span(4,6) ).t(), ones<mat>(1,1) ) ;
		og__.slice(i) = ogh_.subvec(0,2);
	}

	H.clear();}