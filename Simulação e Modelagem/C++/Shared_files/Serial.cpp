#include "Serial.h"

Mecanismo::Mecanismo(uint dof){
	this->dof = dof;
	dy = new Dy(dof);
	q0_.zeros(dof);
	q1_.zeros(dof);
}

Mecanismo::~Mecanismo(){
	delete dy;
	q0_.clear();
	q1_.clear();
}

Dy* Mecanismo::Doit(vec q0_, vec q1_){
	this->q0_ = q0_;
	this->q1_ = q1_;
	return dy;
}

Serial::Serial(int dof, vec l_, vec lg_, vec m_, cube I__, vec g_, mat (*fDH)(vec, vec, vec)):Mecanismo(dof){
	this->l_ = l_;
	this->lg_= lg_;
	this->m_= m_;
	this->I__= I__;
	this->g_ = g_;
	this->fDH = fDH;
	Ig__.zeros(3,3,dof);
	Hr__.zeros(4,4,dof);
	H__.zeros(4,4,dof);
	z__.zeros(3,1,dof+1);
	z__(2,0,0) = 1;
	o__.zeros(3,1,dof+1);
	og__.zeros(3,1,dof);
	Jv__.zeros(3,dof,dof);
	Jw__.zeros(3,dof,dof);
	//Jw2__.zeros(3,dof,dof);
	Jv_n_.zeros(3,dof);
	Jw_n_.zeros(3,dof);
	Mh_.zeros(dof,dof);
	vh_.zeros(dof);
	gh_.zeros(dof);
	v__.zeros(3,1,dof+1);
	V__.zeros(3,dof+1,dof);
	w__.zeros(3,1,dof+1);
	W__.zeros(3,dof+1,dof);
	//w_rel__.zeros(3,1,dof);
	//w_arr__.zeros(3,1,dof);
	//dw_co2__.zeros(3,1,dof);
	//a_cen__.zeros(3,1,dof);
	//a_cor__.zeros(3,1,dof);
	a_co__.zeros(3,1,dof);
	dw_co__.zeros(3,1,dof);
	a_co_n_.zeros(3);
	dw_co_n_.zeros(3);
}

Serial::~Serial(){
	l_.clear() ;
	lg_.clear();
	Ig__.clear();
	Hr__.clear();
	H__.clear();
	z__.clear();
	o__.clear();
	og__.clear();
	Jv__.clear();
	Jw__.clear();
	//Jw2__.clear();
	Jv_n_.clear();
	Jw_n_.clear();
	Mh_.clear();
	vh_.clear();
	gh_.clear();
	//delete dy;
	v__.clear();
	V__.clear();
	w__.clear();
	W__.clear();
	//w_rel__.clear();
	//w_arr__.clear();
	//dw_co2__.clear();
	//a_cen__.clear();
	//a_cor__.clear();
	a_co__.clear();
	dw_co__.clear();
	a_co_n_.clear();
	dw_co_n_.clear();
}

Dy* Serial::Doit(vec q0_, vec q1_){
	this->q0_ = q0_;
	this->q1_ = q1_;
	vec ogr_; ogr_.zeros(4);
	vec ogh_; ogh_.zeros(4);
	mat H = this->fDH(q0_, l_, lg_);

	//CINEMATICA DE POSICAO
	for(uint i = 0; i<dof; i++){
		ogr_ = H( span(i,i), span(4,6) ).t();
		Hr__.slice(i) = H_d(H(i,0),H(i,1),H(i,2),H(i,3));
		H__.slice(i) = i==0 ? Hr__.slice(i): H__.slice(i-1)*Hr__.slice(i);
		z__.slice(i+1) = H__( span(0,2), span(2,2), span(i,i) );
		o__.slice(i+1) = H__( span(0,2), span(3,3), span(i,i) );
		ogh_ = H__.slice(i)*join_vert( ogr_, ones<mat>(1,1) ) ;
		og__.slice(i) = ogh_.subvec(0,2);}

	//CALCULO DOS JACOBIANOS
	for(uint i = 0; i<dof; i++)
		for(uint j = 0; j<=i; j++){
			if(H(j,7)==true){
				Jv__( span(0,2), span(j,j), span(i,i) ) = cross(z__.slice(j), og__.slice(i) - o__.slice(j));
				Jw__( span(0,2), span(j,j), span(i,i) ) = z__.slice(j);}
			else
				Jv__( span(0,2), span(j,j), span(i,i) ) = z__.slice(j);}
		//Jw2__.slice(i) = ((H__.slice(i))(span(0,2),span(0,2))).t()*Jw__.slice(i);}

	for(uint i = 0; i<dof; i++)
		Jv_n_( span(0,2), span(i,i)) = H(i,7)==true ? cross(z__.slice(i), o__.slice(dof) - o__.slice(i)) : z__.slice(i);
	Jw_n_ = Jw__.slice(dof-1);

	//CINEMATICA DE VELOCIDADES
	for(uint i = dof; i>=1; i--){
		v__.slice(i-1) = v__.slice(i) + Jv_n_( span(0,2), span(i-1,i-1))*q0_(i-1);
		w__.slice(i-1) = w__.slice(i) + Jw_n_( span(0,2), span(i-1,i-1))*q0_(i-1);
	}

	for(uint i = 0; i<dof; i++)
		for(int j = i; j>=1; j--){
			V__( span(0,2), span(j-1,j-1), span(i,i) ) = V__( span(0,2), span(j,j), span(i,i) ) + Jv__( span(0,2), span(j-1,j-1), span(i,i) )*q0_(j-1);
			W__( span(0,2), span(j-1,j-1), span(i,i) ) = W__( span(0,2), span(j,j), span(i,i) ) + Jw__( span(0,2), span(j-1,j-1), span(i,i) )*q0_(j-1);
		}

	//CALCULO DAS ACELERACOES COMPLEMENTARES
	a_co_n_.zeros();
	dw_co_n_.zeros();
	for(uint i = 0; i<dof; i++){
		a_co_n_ += q1_(i)*cross( Jw_n_(span(0,2), span(i,i) ), q1_(i)*Jv_n_(span(0,2), span(i,i) ) + 2*v__.slice(i+1) );
		dw_co_n_ += q1_(i)*cross( Jw_n_(span(0,2), span(i,i) ),  w__.slice(i+1) ); }

	a_co__.zeros();
	dw_co__.zeros();
	for(uint i = 0; i<dof; i++)
		for(uint j = 0; j<=i; j++){
			a_co__.slice(j) += q1_(j)*cross( Jw__.slice(i)(span(0,2), span(j,j) ), q1_(j)*Jv__.slice(i)(span(0,2), span(j,j) ) + 2*V__.slice(i)( span(0,2), span(j+1,j+1)  ) );
			dw_co__.slice(j) += q1_(j)*cross( Jw__.slice(i)(span(0,2), span(j,j) ),  W__.slice(i)( span(0,2), span(j+1,j+1)  ) ); }

	//DINAMICA
	Mh_.zeros();
	vh_.zeros();
	gh_.zeros();
	for(uint i = 0; i<dof; i++){
		Ig__.slice(i) = ((H__.slice(i))(span(0,2),span(0,2)))*Ig__.slice(i)*((H__.slice(i))(span(0,2),span(0,2))).t();
		Mh_ += m_(i)*Jv__.slice(i).t()*Jv__.slice(i) + Jw__.slice(i).t()*Ig__.slice(i)*Jw__.slice(i);
		vh_+= Jw__.slice(i).t()*(cross(Jw__.slice(i)*q1_, Ig__.slice(i)*Jw__.slice(i)*q1_) + Ig__.slice(i)*dw_co__.slice(i)) + m_(i)*Jv__.slice(i).t()*a_co__.slice(i);
		gh_+= -m_(i)*Jv__.slice(i).t()*g_; }

	dy->Mh_ = Mh_;
	dy->vh_ = vh_;
	dy->gh_ = gh_;

	H.clear();
	ogh_.clear();
	ogr_.clear();

	return dy;
}