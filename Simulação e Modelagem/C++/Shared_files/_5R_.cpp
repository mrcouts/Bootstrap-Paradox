#include "_5R_.h"

_5R_::_5R_(uint dof, Mecanismo *P, Serial **R_, uint nR_):Mecanismo(dof){
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
    field<vec> q0h__(nR_+1);
    field<vec> q1h__(nR_+1);
    q0h__(0) = q0_(span(0,dof-1));
    q1h__(0) = q1_(span(0,dof-1));
    uint soma = dof;
    for(uint i=0; i<nR_; i++){
        soma += R_[i]->dof;
        q0h__(i+1) = q0_(span(soma-R_[i]->dof,soma-1));
        q1h__(i+1) = q1_(span(soma-R_[i]->dof,soma-1));
    }

    P->Doit(q0h__(0), q1h__(0));
    for(uint i=0; i<nR_; i++) R_[i]->Doit(q0h__(i+1), q1h__(i+1));
    //M_, v_ e g_
    M_ = join_diag(P->dy->Mh_, join_diag(R_[0]->Mh_, R_[1]->Mh_) );
    v_ = join_vert(P->dy->vh_, join_vert(R_[0]->vh_, R_[1]->vh_) );
    g_ = join_vert(P->dy->gh_, join_vert(R_[0]->gh_, R_[1]->gh_) );
    // _q_, Ah_, Ao_, A_, b_, C_ e Z_
    mat A = join_diag( Roty(0)(span(0,1),span(0,2)), Roty(PI)(span(0,1),span(0,2)) );
    vec b = {0.05,0,-0.05,0};
    _q_ = join_vert(q0h__(0), q0h__(0)) - A*join_vert(R_[0]->o__.slice(R_[0]->dof), R_[1]->o__.slice(R_[1]->dof)) - b;
    Ah_ = join_vert(eye(dof,dof),eye(dof,dof));
    Ao_ = -A*join_diag(R_[0]->Jv_n_,R_[1]->Jv_n_);
    A_ = join_horiz(Ah_, Ao_);
    b_ = A*join_vert(R_[0]->a_co_n_, R_[1]->a_co_n_);
    C_ = join_vert( eye(dof,dof), -solve(Ao_,Ah_) );
    Z_ = join_vert(C_.row(2), C_.row(4));
    // Mh_, vh_ e gh_
    dy->Mh_ = solve(Z_.t(), C_.t()*M_*C_);
    dy->vh_ = solve(Z_.t(), C_.t()*( M_*join_vert(zeros(dof), solve(Ao_, b_) ) + v_ ));
    dy->gh_ = solve(Z_.t(), C_.t()*g_);

    A.clear();
    b.clear();
    q0h__.clear();
    q1h__.clear();

    return dy;
}