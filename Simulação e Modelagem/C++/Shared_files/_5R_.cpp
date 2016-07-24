#include "_5R_.h"

_5R_::_5R_():Mecanismo(2){
    dof = 2;
    nq = 6;
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
    Mh_.zeros(dof,dof);
    vh_.zeros(dof);
    gh_.zeros(dof);
    H_.zeros(dof,dof);
    h_.zeros(dof);

    cube I__; I__.zeros(3,3,2);
    I__.slice(0) << 0 << 0        << 0        << endr
                 << 0 << 171.6e-6 << 0        << endr
                 << 0 << 0        << 171.6e-6 << endr;

    I__.slice(1) << 0 << 0        << 0        << endr
                 << 0 << 320.6e-6 << 0        << endr
                 << 0 << 0        << 320.6e-6 << endr;

    RR1 = new Serial(2, {0.12, 0.15}, {0.5*0.12, 0.5*0.15},{0.143, 0.171}, I__ , {0, -9.8,0}, &fDH_RR);
    RR2 = new Serial(2, {0.12, 0.15}, {0.5*0.12, 0.5*0.15},{0.143, 0.171}, I__ , {0, -9.8,0}, &fDH_RR);
    P = new Mecanismo(dof); 
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
    RR1->Doit(q0h1_, q1h1_);
    RR2->Doit(q0h2_, q1h2_);
    //M_, v_ e g_
    M_ = join_diag(P->dy->Mh_, join_diag(RR1->Mh_, RR2->Mh_) );
    v_ = join_vert(P->dy->vh_, join_vert(RR1->vh_, RR2->vh_) );
    g_ = join_vert(P->dy->gh_, join_vert(RR1->gh_, RR2->gh_) );
    // _q_, Ah_, Ao_, A_, b_, C_ e Z_
    mat A = join_diag( Roty(0)(span(0,1),span(0,2)), Roty(PI)(span(0,1),span(0,2)) );
    vec b = {0.05,0,-0.05,0};
    _q_ = join_vert(q0h_, q0h_) - A*join_vert(RR1->o__.slice(RR1->dof), RR2->o__.slice(RR2->dof)) - b;
    Ah_ = join_vert(eye(dof,dof),eye(dof,dof));
    Ao_ = -A*join_diag(RR1->Jv_n_,RR2->Jv_n_);
    A_ = join_horiz(Ah_, Ao_);
    b_ = A*join_vert(RR1->a_co_n_, RR2->a_co_n_);
    C_ = join_vert( eye(dof,dof), -solve(Ao_,Ah_) );
    Z_ = join_vert(C_.row(2), C_.row(4));
    // Mh_, vh_ e gh_
    Mh_ = C_.t()*M_*C_;
    vh_ = C_.t()*( M_*join_vert(zeros(dof), solve(Ao_, b_) ) + v_ );
    gh_ = C_.t()*g_;
    dy->Mh_ = solve(Z_.t(),Mh_);
    dy->vh_ = solve(Z_.t(),vh_);
    dy->gh_ = solve(Z_.t(),gh_);
    // H_ e h_
    //H_ = solve(Z_.t(), Mh_);
    //h_ = solve(Z_.t(), vh_ + gh_);
    return dy;
}