#include "Serial.h"

Mecanismo::Mecanismo(int dof){
	dy = new Dy(dof); }

Mecanismo::~Mecanismo(){
	delete dy; }

Dy* Mecanismo::Doit(vec q0_, vec q1_){
	return dy;}

Serial::Serial(int dof, vec l_, vec lg_, vec m_, cube I__, vec g_, mat (*fDH)(vec, vec, vec)):Mecanismo(dof){
	this->dof = dof;
	this->l_ = l_;
	this->lg_= lg_;
	this->m_= m_;
	this->I__= I__;
	this->g_ = g_;
	this->fDH = fDH;
	Hr__.zeros(4,4,dof);
	H__.zeros(4,4,dof);
	z__.zeros(3,1,dof+1);
	z__(2,0,0) = 1;
	o__.zeros(3,1,dof+1);
	og__.zeros(3,1,dof);
	Jv__.zeros(3,dof,dof);
	Jw__.zeros(3,dof,dof);
	Jw2__.zeros(3,dof,dof);
	Jv_n_.zeros(3,dof);
	Jw_n_.zeros(3,dof);
	Mh_.zeros(dof,dof);
	vh_.zeros(dof);
	gh_.zeros(dof);
	//dy = new Dy(dof);
	w_rel__.zeros(3,1,dof);
	w_arr__.zeros(3,1,dof);
	dw_co__.zeros(3,1,dof);
	dw_co2__.zeros(3,1,dof);
	a_cen__.zeros(3,1,dof);
	a_cor__.zeros(3,1,dof);
	a_co__.zeros(3,1,dof);
	dw_co_n_.zeros(3);
	a_co_n_.zeros(3); }

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
	Jw2__.clear();
	Jv_n_.clear();
	Jw_n_.clear();
	Mh_.clear();
	vh_.clear();
	gh_.clear();
	//delete dy;
	w_rel__.clear();
	w_arr__.clear();
	dw_co__.clear();
	dw_co2__.clear();
	a_cen__.clear();
	a_cor__.clear();
	a_co__.clear();
	dw_co_n_.clear();
	a_co_n_.clear(); }

Dy* Serial::Doit(vec q0_, vec q1_){
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
				Jv__( span(0,2), span(j,j), span(i,i) ) = z__.slice(j);}
		Jw2__.slice(i) = ((H__.slice(i))(span(0,2),span(0,2))).t()*Jw__.slice(i);}

	for(int i = 0; i<dof; i++)
		Jv_n_( span(0,2), span(i,i)) = H(i,7)==true ? cross(z__.slice(i), o__.slice(dof) - o__.slice(i)) : z__.slice(i);
	Jw_n_ = Jw__.slice(dof-1);

	Mh_.zeros();
	vh_.zeros();
	gh_.zeros();
	for(int i = 0; i<dof; i++){
		Mh_ += m_(i)*Jv__.slice(i).t()*Jv__.slice(i) + Jw2__.slice(i).t()*I__.slice(i)*Jw2__.slice(i);
		vh_+= Jw2__.slice(i).t()*cross(Jw2__.slice(i)*q1_, I__.slice(i)*Jw2__.slice(i)*q1_);
		gh_+= -m_(i)*Jv__.slice(i).t()*g_; }

	for(int i = 0; i<dof; i++)
		w_rel__.slice(i) = Jw_n_( span(0,2), span(i,i) )*q1_(i);

	for(int i = 1; i<dof; i++) {
		w_arr__.slice(i) = w_arr__.slice(i-1) + w_rel__.slice(i-1);
		dw_co__.slice(i) = dw_co__.slice(i-1) + cross(w_arr__.slice(i), w_rel__.slice(i));
		dw_co2__.slice(i) = ((H__.slice(i))(span(0,2),span(0,2))).t()*dw_co__.slice(i); }

	a_cen__.zeros();
	a_cor__.zeros();
	vec v_rel; v_rel.zeros(3);
	for(int i = 0; i<dof; i++){
		for(int j = 0; j<=i; j++){
			v_rel = Jv__(span(0,2),span(j,j),span(i,i))*q1_(j);
			a_cor__.slice(i) += 2*cross(w_arr__.slice(j), v_rel );
			if(H(j,7)==true) 
				a_cen__.slice(i) += cross(w_rel__.slice(j), cross(w_rel__.slice(j), og__.slice(i) - o__.slice(j)) ); }}
	a_co__ = a_cen__ + a_cor__;

	for(int i = 0; i<dof; i++)
		vh_ += m_(i)*Jv__.slice(i).t()*a_co__.slice(i) + Jw2__.slice(i).t()*I__.slice(i)*dw_co2__.slice(i);
	dy->Mh_ = Mh_;
	dy->vh_ = vh_;
	dy->gh_ = gh_;

	dw_co_n_ = dw_co__.slice(dof-1);

	a_co_n_.zeros();
	vec a_cor_n_; a_cor_n_.zeros(3); 
	vec a_cen_n_; a_cen_n_.zeros(3); 
	v_rel.zeros(3);
	for(int j = 0; j<dof; j++){
		v_rel = Jv_n_(span(0,2),span(j,j))*q1_(j);
		a_cor_n_ += 2*cross(w_arr__.slice(j), v_rel );
		if(H(j,7)==true) 
			a_cen_n_ += cross(w_rel__.slice(j), cross(w_rel__.slice(j), o__.slice(dof) - o__.slice(j)) ); }
	a_co_n_ = a_cen_n_ + a_cor_n_;


	H.clear();
	ogh_.clear();
	a_cor_n_.clear();
	a_cen_n_.clear();
	v_rel.clear();

	return dy;}