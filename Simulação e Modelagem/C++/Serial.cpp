#include "Serial.h"

Serial::Serial(int dof, vec l_, vec lg_, mat (*fDH)(vec, vec, vec)){
	this->dof = dof;
	this->l_ = l_;
	this->lg_= lg_;
	this->fDH = fDH;
	Hr__.zeros(4,4,dof);
	H__.zeros(4,4,dof);
	z__.zeros(3,1,dof+1);
	z__(2,0,0) = 1;
	o__.zeros(3,1,dof+1);
	og__.zeros(3,1,dof);
	Jv__.zeros(3,dof,dof);
	Jw__.zeros(3,dof,dof);
	Jv_n_.zeros(3,dof);
	Jw_n_.zeros(3,dof); }

Serial::~Serial(){
	l_.clear() ;
	lg_.clear();
	Hr__.clear();
	H__.clear();
	z__.clear();
	o__.clear();
	og__.clear();
	Jv__.clear();
	Jw__.clear();
	Jv_n_.clear();
	Jw_n_.clear(); }

void Serial::Doit(vec q0_){
	vec ogh_; ogh_.zeros(4);
	mat H = this->fDH(q0_, l_, lg_);
	for(int i = 0; i<dof; i++){
		Hr__.slice(i) = H_d(H(i,0),H(i,1),H(i,2),H(i,3));
		H__.slice(i) = i==0 ? Hr__.slice(i): H__.slice(i-1)*Hr__.slice(i);
		z__.slice(i+1) = H__( span(0,2), span(2,2), span(i,i) );
		o__.slice(i+1) = H__( span(0,2), span(3,3), span(i,i) );
		ogh_ = H__.slice(i)*join_vert( H( span(i,i), span(4,6) ).t(), ones<mat>(1,1) ) ;
		og__.slice(i) = ogh_.subvec(0,2);}


	for(int i = 0; i<dof; i++){
		for(int j = 0; j<=i; j++){
			if(H(j,7)==true){
				Jv__( span(0,2), span(j,j), span(i,i) ) = cross(z__.slice(j), og__.slice(i) - o__.slice(j));
				Jw__( span(0,2), span(j,j), span(i,i) ) = z__.slice(j);}
			else
				Jv__( span(0,2), span(j,j), span(i,i) ) = z__.slice(j);}}

	for(int i = 0; i<dof; i++)
		Jv_n_( span(0,2), span(i,i)) = H(i,7)==true ? cross(z__.slice(i), o__.slice(dof) - o__.slice(i)) : z__.slice(i);

	Jw_n_ = Jw__.slice(dof-1);

	H.clear();
	ogh_.clear();}