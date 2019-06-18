#include "Parallel.h"

Parallel::Parallel(uint dof, Mecanismo *P, Serial **R_, uint nR_, Col<uint> u_nzi_, mat D_, mat E_, mat F_, vec f_):Mecanismo(dof){
    this->Construct(dof, P, R_, nR_, u_nzi_, D_, E_, F_, f_);
    Ah_ = D_;
    caso = 1;
}

Parallel::Parallel(uint dof, Mecanismo *P, Serial **R_, uint nR_, Col<uint> u_nzi_, mat D_, mat E_, mat F_, vec f_, mat P_, mat Q_, mat S_):Mecanismo(dof){
    this->Construct(dof, P, R_, nR_, u_nzi_, D_, E_, F_, f_);
    this->P_ = P_;
    this->Q_ = Q_;
    this->S_ = S_;
    Ah_ = join_vert(D_, P_);
    caso = 2;
}

void Parallel::Construct(uint dof, Mecanismo *P, Serial **R_, uint nR_, Col<uint> u_nzi_, mat D_, mat E_, mat F_, vec f_){
	this->P = P;
    this->R_ = R_;
    this->nR_= nR_;
    this->u_nzi_ = u_nzi_;
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

    q0__      = new vec* [nR_]; for(uint i = 0; i<nR_; i++)      q0__[i] = &(R_[i]->q0_);
    q1__      = new vec* [nR_]; for(uint i = 0; i<nR_; i++)      q1__[i] = &(R_[i]->q1_);
    o__       = new vec* [nR_]; for(uint i = 0; i<nR_; i++)       o__[i] = &(R_[i]->x_);
    Jv_n__    = new mat* [nR_]; for(uint i = 0; i<nR_; i++)    Jv_n__[i] = &(R_[i]->Jv_n_);
    Jw_n__    = new mat* [nR_]; for(uint i = 0; i<nR_; i++)    Jw_n__[i] = &(R_[i]->Jw_n_);
    a_co_n__  = new vec* [nR_]; for(uint i = 0; i<nR_; i++)  a_co_n__[i] = &(R_[i]->a_co_n_);
    dw_co_n__ = new vec* [nR_]; for(uint i = 0; i<nR_; i++) dw_co_n__[i] = &(R_[i]->dw_co_n_);

    this->D_ = D_;
    this->E_ = E_;
    this->F_ = F_;
    this->f_ = f_;
}

Parallel::~Parallel(){
    M_.clear();
    v_.clear();
    g_.clear();
    _q_.clear();
    Ah_.clear();
    Ao1_.clear();
    Ao2_.clear();
    Ao_.clear();
    A_.clear();
    b1_.clear();
    b2_.clear();
    b_.clear();
    C_.clear();
    Z_.clear();
}

void Parallel::DoitKin(vec q0_, vec q1_){
    this->q0_ = q0_;
    this->q1_ = q1_;

    P->Doit(q0_(span(0,P->dof-1)), q1_(span(0,P->dof-1)));
    uint soma = P->dof;
    for(uint i=0; i<nR_; i++){
        soma += R_[i]->dof;
        R_[i]->Doit(q0_(span(soma-R_[i]->dof,soma-1)), q1_(span(soma-R_[i]->dof,soma-1)) );
    }

    // _q_, Ah_, Ao_, A_, b_, C_ e Z_
    _q_ = D_*P->q0_ - E_*join_vert(o__, nR_) - F_*join_vert(q0__, nR_) - f_;
    Ao1_ = -E_*join_diag(Jv_n__, nR_) - F_;
    b1_ = E_*join_vert(a_co_n__, nR_);
    switch(caso){
        case 1:
            Ao_ = Ao1_;            
            b_ = b1_;            
            break;
        case 2:
            Ao2_ = -Q_*join_diag(Jw_n__, nR_) - S_;
            b2_ = Q_*join_vert(dw_co_n__, nR_);
            Ao_ = join_vert(Ao1_, Ao2_);
            b_ = join_vert(b1_, b2_);
            break;
    }
    A_ = join_horiz(Ah_, Ao_);
}

Dy* Parallel::Doit(vec q0_, vec q1_){
    this->DoitKin(q0_, q1_);

    C_ = join_vert( eye(dof,dof), -solve(Ao_,Ah_) );
    mat Aux_ = C_.row(u_nzi_(0)); for(uint i=1; i<u_nzi_.n_rows; i++) Aux_ = join_vert(Aux_, C_.row(u_nzi_(i)));
    Z_ = Aux_;

    //M_, v_ e g_
    M_ = join_diag(M__, nR_+1);
    v_ = join_vert(v__, nR_+1);
    g_ = join_vert(g__, nR_+1);
    // Mh_, vh_ e gh_
    dy->Mh_ = solve(Z_.t(), C_.t()*M_*C_);
    dy->vh_ = solve(Z_.t(), C_.t()*( M_*join_vert(zeros(dof), solve(Ao_, b_) ) + v_ ));
    dy->gh_ = solve(Z_.t(), C_.t()*g_);
    //mat invZ_ = inv(Z_);
    //dy->Mh_ = invZ_.t()*(C_.t()*M_*C_)*invZ_;
    //dy->vh_ = invZ_.t()*(C_.t()*( M_*join_vert(zeros(dof), solve(Ao_, b_) ) + v_ ));
    //dy->gh_ = invZ_.t()*(C_.t()*g_);
    return dy;
}