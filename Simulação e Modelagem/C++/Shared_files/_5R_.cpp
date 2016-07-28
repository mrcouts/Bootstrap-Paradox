#include "_5R_.h"

_5R_::_5R_(uint dof, Mecanismo *P, Serial **R_, uint nR_, mat D_, mat E_, mat F_, vec f_):Mecanismo(dof){
    this->P = P;
    this->R_ = R_;
    this->nR_= nR_;
    nq = P->dof; for(uint i = 0; i<nR_; i++) nq += R_[i]->dof;
    M_.zeros(nq,nq);
    v_.zeros(nq);
    g_.zeros(nq);
    _q_.zeros(nq-dof);
    Ah_.zeros(nq-dof,dof);
    Ao_.zeros(nq-dof,nq-dof);
    A_.zeros(nq-dof,nq);
    b_.zeros(nq-dof);
    C_.zeros(nq,dof);
    Z_.zeros(dof,dof);

    M__ = new mat* [nR_+1];
    M__[0] = &(P->dy->Mh_); for(uint i = 0; i<nR_; i++) M__[i+1] = &(R_[i]->Mh_);
    v__ = new vec* [nR_+1];
    v__[0] = &(P->dy->vh_); for(uint i = 0; i<nR_; i++) v__[i+1] = &(R_[i]->vh_);
    g__ = new vec* [nR_+1];
    g__[0] = &(P->dy->gh_); for(uint i = 0; i<nR_; i++) g__[i+1] = &(R_[i]->gh_);

    q0__     = new vec* [nR_]; for(uint i = 0; i<nR_; i++)     q0__[i] = &(R_[i]->q0_);
    q1__     = new vec* [nR_]; for(uint i = 0; i<nR_; i++)     q1__[i] = &(R_[i]->q1_);
    o__      = new mat* [nR_]; for(uint i = 0; i<nR_; i++)      o__[i] = &(R_[i]->o__.slice(R_[i]->dof));
    Jv_n__   = new mat* [nR_]; for(uint i = 0; i<nR_; i++)   Jv_n__[i] = &(R_[i]->Jv_n_);
    a_co_n__ = new mat* [nR_]; for(uint i = 0; i<nR_; i++) a_co_n__[i] = &(R_[i]->a_co_n_);

    this->D_ = D_;
    this->E_ = E_;
    this->F_ = F_;
    this->f_ = f_;

    Ah_ = D_;
}

_5R_::~_5R_(){
    M_.clear();
    v_.clear();
    g_.clear();
    _q_.clear();
    Ah_.clear();
    Ao_.clear();
    A_.clear();
    b_.clear();
    C_.clear();
    Z_.clear();
}

Dy* _5R_::Doit(vec q0_, vec q1_){
	this->q0_ = q0_;
	this->q1_ = q1_;

    P->Doit(q0_(span(0,P->dof-1)), q1_(span(0,P->dof-1)));
    uint soma = P->dof;
    for(uint i=0; i<nR_; i++){
        soma += R_[i]->dof;
        R_[i]->Doit(q0_(span(soma-R_[i]->dof,soma-1)), q1_(span(soma-R_[i]->dof,soma-1)) );
    }

    //M_, v_ e g_
    M_ = join_diag(M__, nR_+1);
    v_ = join_vert(v__, nR_+1);
    g_ = join_vert(g__, nR_+1);
    // _q_, Ah_, Ao_, A_, b_, C_ e Z_
    _q_ = D_*P->q0_ - E_*join_vert(o__, nR_) - F_*join_vert(q0__, nR_) - f_;
    Ao_ = -E_*join_diag(Jv_n__, nR_) - F_;
    A_ = join_horiz(Ah_, Ao_);
    b_ = E_*join_vert(a_co_n__, nR_);
    C_ = join_vert( eye(dof,dof), -solve(Ao_,Ah_) );
    Z_ = join_vert(C_.row(2), C_.row(4));
    // Mh_, vh_ e gh_
    dy->Mh_ = solve(Z_.t(), C_.t()*M_*C_);
    dy->vh_ = solve(Z_.t(), C_.t()*( M_*join_vert(zeros(dof), solve(Ao_, b_) ) + v_ ));
    dy->gh_ = solve(Z_.t(), C_.t()*g_);

    return dy;
}