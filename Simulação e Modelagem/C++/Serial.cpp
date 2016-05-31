#include "Serial.h"

Serial::Serial(int dof, vec l_, vec lg_, mat (*fDH)(vec, vec, vec)){
	this->dof = dof;
	this->l_ = l_;
	this->lg_= lg_;
	this->Hr__.zeros(4,4,dof);
	this->H__.zeros(4,4,dof);
	this->fDH = fDH; }

Serial::~Serial(){
	l_.clear() ;
	lg_.clear();
	Hr__.clear();
	H__.clear(); }

void Serial::Doit(vec q0_){
	mat H;
	for(int i = 0; i<dof; i++){
		H = this->fDH(q0_, l_, lg_);
		Hr__.slice(i) = H_d(H(i,0),H(i,1),H(i,2),H(i,3));
		if(i == 0)
			H__.slice(i) = Hr__.slice(i);
		else
			H__.slice(i) = H__.slice(i-1)*Hr__.slice(i);}
	H.clear();}