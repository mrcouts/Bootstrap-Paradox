#include "_5R_.h"

_5R_::_5R_(uint dof, Mecanismo *P, Serial **RR_):Mecanismo(dof){
    nq = P->dof + RR_[0]->dof + RR_[1]->dof;
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
    this->P = P;
    this->RR_ = RR_;
}

_5R_::~_5R_(){}

Dy* _5R_::Doit(vec q0_, vec q1_){
   	vec q0h_ = q0_(span(0,1));
   	vec q1h_ = q1_(span(0,1));
   	vec q0h1_ = q0_(span(2,3));
   	vec q1h1_ = q1_(span(2,3));
   	vec q0h2_ = q0_(span(4,5));
   	vec q1h2_ = q1_(span(4,5));

    P->Doit(q0h_, q1h_);
    RR_[0]->Doit(q0h1_, q1h1_);
    RR_[1]->Doit(q0h2_, q1h2_);
    //M_, v_ e g_
    M_ = join_diag(P->dy->Mh_, join_diag(RR_[0]->Mh_, RR_[1]->Mh_) );
    v_ = join_vert(P->dy->vh_, join_vert(RR_[0]->vh_, RR_[1]->vh_) );
    g_ = join_vert(P->dy->gh_, join_vert(RR_[0]->gh_, RR_[1]->gh_) );
    // _q_, Ah_, Ao_, A_, b_, C_ e Z_
    mat A = join_diag( Roty(0)(span(0,1),span(0,2)), Roty(PI)(span(0,1),span(0,2)) );
    vec b = {0.05,0,-0.05,0};
    _q_ = join_vert(q0h_, q0h_) - A*join_vert(RR_[0]->o__.slice(RR_[0]->dof), RR_[1]->o__.slice(RR_[1]->dof)) - b;
    Ah_ = join_vert(eye(dof,dof),eye(dof,dof));
    Ao_ = -A*join_diag(RR_[0]->Jv_n_,RR_[1]->Jv_n_);
    A_ = join_horiz(Ah_, Ao_);
    b_ = A*join_vert(RR_[0]->a_co_n_, RR_[1]->a_co_n_);
    C_ = join_vert( eye(dof,dof), -solve(Ao_,Ah_) );
    Z_ = join_vert(C_.row(2), C_.row(4));
    // Mh_, vh_ e gh_
    dy->Mh_ = solve(Z_.t(), C_.t()*M_*C_);
    dy->vh_ = solve(Z_.t(), C_.t()*( M_*join_vert(zeros(dof), solve(Ao_, b_) ) + v_ ));
    dy->gh_ = solve(Z_.t(), C_.t()*g_);
    return dy;
}