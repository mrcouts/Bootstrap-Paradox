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

CorpoRigido::CorpoRigido(double m, mat33 I_, vec3 g_){
	dof = 6;
	dy = new Dy(6);
	this->m = m;
	this->I_ = I_;
	this->g_ = g_;
}

CorpoRigido::~CorpoRigido(){
	delete dy;
}

Dy* CorpoRigido::Doit(vec7 q_, vec6 p_){
	x_ = q_.subvec(0,2);
	qu_ = q_.subvec(3,6);
	qu_ = qu_/norm(qu_);
	this->q_ = join_vert(x_,qu_);

	v_ = p_.subvec(0,2);
	w_ = p_.subvec(3,5);
	this->p_ = p_;

	Ce_ << -qu_(1) << -qu_(2) << -qu_(3) << endr
       <<   qu_(0) <<  qu_(3) << -qu_(2) << endr
       <<  -qu_(3) <<  qu_(0) <<  qu_(1) << endr
       <<   qu_(2) << -qu_(1) <<  qu_(0) << endr;

    dq_p_ = join_diag( eye(3,3), 0.5*Ce_ );
    p_dq_ = join_diag( eye(3,3), 2*Ce_.t() );

    R_ << 1 - 2*(qu_(2)*qu_(2) + qu_(3)*qu_(3))  <<     2*(qu_(1)*qu_(2) - qu_(3)*qu_(0)) <<     2*(qu_(1)*qu_(3) + qu_(2)*qu_(0)) << endr
       <<     2*(qu_(1)*qu_(2) + qu_(3)*qu_(0))  << 1 - 2*(qu_(1)*qu_(1) + qu_(3)*qu_(3)) <<     2*(qu_(2)*qu_(3) - qu_(1)*qu_(0)) << endr
       <<     2*(qu_(1)*qu_(3) - qu_(2)*qu_(0))  <<     2*(qu_(1)*qu_(0) + qu_(2)*qu_(3)) << 1 - 2*(qu_(1)*qu_(1) + qu_(2)*qu_(2)) << endr;

    Ig_ = R_*I_*R_.t();
   	dy->Mh_ = join_diag( m*eye(3,3), Ig_ );
   	dy->vh_ = join_vert( zeros(3,1), cross(w_, Ig_*w_) );
   	dy->gh_ = join_vert(-m*g_, zeros(3,1));

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
	z__.zeros(3,dof+1);
	z__(2,0) = 1;
	o__.zeros(3,dof+1);
	x_.zeros(3);
	og__.zeros(3,dof);
	Jv__.zeros(3,dof,dof);
	Jw__.zeros(3,dof,dof);
	Jv_n_.zeros(3,dof);
	Jw_n_.zeros(3,dof);
	Mh_.zeros(dof,dof);
	vh_.zeros(dof);
	gh_.zeros(dof);
	v__.zeros(3,dof+1);
	v___.zeros(3,dof+1,dof);
	w__.zeros(3,dof+1);
	a_co_n_i_.zeros(3,dof);
	dw_co_n_i_.zeros(3,dof);
	a_co_n_.zeros(3);
	dw_co_n_.zeros(3);
	a_co_ij__.zeros(3,dof,dof);
	a_co__.zeros(3,dof);
	dw_co__.zeros(3,dof);
}

Serial::~Serial(){
	l_.clear() ;
	lg_.clear();
	Ig__.clear();
	Hr__.clear();
	H__.clear();
	z__.clear();
	o__.clear();
	x_.clear();
	og__.clear();
	Jv__.clear();
	Jw__.clear();
	Jv_n_.clear();
	Jw_n_.clear();
	Mh_.clear();
	vh_.clear();
	gh_.clear();
	//delete dy;
	v__.clear();
	v___.clear();
	w__.clear();
	a_co_n_.clear();
	dw_co_n_.clear();
	a_co_n_i_.clear();
	dw_co_n_i_.clear();
	a_co__.clear();
	dw_co__.clear();
	a_co_ij__.clear();
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
		z__.col(i+1) = H__.slice(i).col(2).rows(0,2);
		o__.col(i+1) = H__.slice(i).col(3).rows(0,2);
		ogh_ = H__.slice(i)*join_vert( ogr_, ones<mat>(1,1) ) ;
		og__.col(i) = ogh_.rows(0,2);}
	x_ = o__.col(dof);

	//CALCULO DOS JACOBIANOS
	Jw_n_.zeros();
	for(uint i = 0; i<dof; i++){
		if(H(i,7)==true){
			Jv_n_.col(i) = cross(z__.col(i), o__.col(dof) - o__.col(i));
			Jw_n_.col(i) = z__.col(i);}
		else
			Jv_n_.col(i) = z__.col(i);
		Jw__.slice(i) = Jw_n_;
		for(uint j = 0; j<=i; j++){
			if(H(j,7)==true)
				Jv__.slice(i).col(j) = cross(z__.col(j), og__.col(i) - o__.col(j));
			else
				Jv__.slice(i).col(j) = z__.col(j);}}
		
	//CINEMATICA DE VELOCIDADES
	for(uint i = dof; i>=1; i--){
		v__.col(i-1) = v__.col(i) + Jv_n_.col(i-1)*q0_(i-1);
		w__.col(i-1) = w__.col(i) + Jw_n_.col(i-1)*q0_(i-1);}

	for(uint i = 0; i<dof; i++)
		for(int j = i; j>=1; j--)
			v___.slice(i).col(j-1) = v___.slice(i).col(j) + Jv__.slice(i).col(j-1)*q0_(j-1);

	//CALCULO DAS ACELERACOES COMPLEMENTARES
	a_co_n_.zeros();
	dw_co_n_.zeros();
	a_co__.zeros();
	for(uint i = 0; i<dof; i++){
		a_co_n_i_.col(i)  = q1_(i)*cross( Jw_n_.col(i),  q1_(i)*Jv_n_.col(i) + 2*v__.col(i+1) );
		dw_co_n_i_.col(i) = q1_(i)*cross( Jw_n_.col(i),  w__.col(i+1) );
		a_co_n_  += a_co_n_i_.col(i);
		dw_co_n_ += dw_co_n_i_.col(i);
		dw_co__.col(i) = dw_co_n_;
		for(uint j = 0; j<=i; j++){
			a_co_ij__.slice(i).col(j) = q1_(j)*cross( Jw_n_.col(j), q1_(j)*Jv__.slice(i).col(j) + 2*v___.slice(i).col(j+1) );
			a_co__.col(j)  += a_co_ij__.slice(i).col(j); }}

	//DINAMICA
	Mh_.zeros();
	vh_.zeros();
	gh_.zeros();
	for(uint i = 0; i<dof; i++){
		Ig__.slice(i) = ((H__.slice(i))(span(0,2),span(0,2)))*Ig__.slice(i)*((H__.slice(i))(span(0,2),span(0,2))).t();
		Mh_ += m_(i)*Jv__.slice(i).t()*Jv__.slice(i) + Jw__.slice(i).t()*Ig__.slice(i)*Jw__.slice(i);
		vh_+= Jw__.slice(i).t()*(cross(Jw__.slice(i)*q1_, Ig__.slice(i)*Jw__.slice(i)*q1_) + Ig__.slice(i)*dw_co__.col(i)) + m_(i)*Jv__.slice(i).t()*a_co__.col(i);
		gh_+= -m_(i)*Jv__.slice(i).t()*g_; }

	dy->Mh_ = Mh_;
	dy->vh_ = vh_;
	dy->gh_ = gh_;

	H.clear();
	ogh_.clear();
	ogr_.clear();

	return dy;
}