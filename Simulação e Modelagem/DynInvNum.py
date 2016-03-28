from InKinNum import *
from DynamicsDeduction import *

PARAMETROS = {symbols('l_0'): par[l0],
              symbols('l_1'): par[l1],
              symbols('l_2'): par[l2],
              symbols('lg_1'): 0.5*par[l1],
              symbols('lg_2'): 0.5*par[l2],
              symbols('m_1'): 0.143,
              symbols('m_2'): 0.171,
              symbols('Jz_1'): (0.143*par[l1]**2)/12.0,
              symbols('Jz_2'): (0.171*par[l2]**2)/12.0,
              symbols('g'): 9.8         
}
                  
rep = lambda vec1,vec2: [(vec1[i],vec2[i]) for i in xrange(len(vec1))]

Ah_s = Ah_.subs(PARAMETROS)
Ao_s = Ao_.subs(PARAMETROS)
M_s = M_.subs(PARAMETROS)
g_s = g_.subs(PARAMETROS)
v_s = v_.subs(PARAMETROS)
b_s = b_.subs(PARAMETROS)
#C0_s = P.C_.subs(PARAMETROS)
#C1_s = RR1.C_.subs(PARAMETROS)
#C2_s = RR2.C_.subs(PARAMETROS)
#Jh_s = Jh_.subs(PARAMETROS)
#Jo_s = Jo_.subs(PARAMETROS)
#bch_s = bch_.subs(PARAMETROS)
#b1_s = RR1.b_.subs(PARAMETROS)
#b2_s = RR2.b_.subs(PARAMETROS)

Ah_n= lambda q_n:       SMatrix( (Ah_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), Ah_.rowl_, Ah_.coll_)
Ao_n= lambda q_n:       SMatrix( (Ao_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), Ao_.rowl_, Ao_.coll_)
M_n = lambda q_n:       SMatrix( ( M_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), M_.rowl_, M_.coll_)
g_n = lambda q_n:       SMatrix( ( g_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), g_.rowl_, g_.coll_)
v_n = lambda q_n, p_n:  SMatrix( ( v_s.M_.subs(rep( list(p_.M_), p_n)).subs(rep( list(q_.M_) ,q_n))).evalf(),  v_.rowl_,  v_.coll_)
b_n = lambda q_n, p_n:  SMatrix( ( b_s.M_.subs(rep( list(p_.M_), p_n)).subs(rep( list(q_.M_) ,q_n))).evalf(),  b_.rowl_,  b_.coll_)
#C0_n =lambda q_n:       SMatrix( (C0_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), C0_s.rowl_, C0_s.coll_)
#C1_n =lambda q_n:       SMatrix( (C1_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), C1_s.rowl_, C1_s.coll_)
#C2_n =lambda q_n:       SMatrix( (C2_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), C2_s.rowl_, C2_s.coll_)
#Jh_n =lambda q_n:       SMatrix( (Jh_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), Jh_s.rowl_, Jh_s.coll_)
#Jo_n =lambda q_n:       SMatrix( (Jo_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), Jo_s.rowl_, Jo_s.coll_)
#bch_n=lambda q_n, p_n:  SMatrix((bch_s.M_.subs(rep( list(p_.M_), p_n)).subs(rep( list(q_.M_) ,q_n))).evalf(),  bch_s.rowl_,  bch_s.coll_)
#b1_n= lambda q_n, p_n:  SMatrix(( b1_s.M_.subs(rep( list(p_.M_), p_n)).subs(rep( list(q_.M_) ,q_n))).evalf(),   b1_s.rowl_,   b1_s.coll_)
#b2_n= lambda q_n, p_n:  SMatrix(( b2_s.M_.subs(rep( list(p_.M_), p_n)).subs(rep( list(q_.M_) ,q_n))).evalf(),   b2_s.rowl_,   b2_s.coll_)
   
def tau_n(qh_t, ph_t, dph_t):
    ph_t = SMatrix( ph_t, ph_.rowl_)
    dph_t = SMatrix( dph_t, ph_.rowl_)
    q_t = fInvKin(qh_t)
    Ah_t = Ah_n(q_t)
    Ao_t = Ao_n(q_t)
    iAo_t = Ao_t.inv()
    C_t = SMatrix(1, ph_.rowl_, ph_.rowl_) - iAo_t*Ah_t
    p_t = C_t*ph_t
    b_t = b_n(q_t, p_t.M_)
    dpo_t = iAo_t*(-1*Ah_t*dph_t + b_t)
    dp_t = dph_t + dpo_t
    Z_t = C_t.T()*U_
    M_t = M_n(q_t)
    v_t = v_n(q_t, p_t.M_)
    g_t = g_n(q_t)    
    tau_t = Z_t.LDLsolve( C_t.T()*( M_t*dp_t + v_t + g_t ) )
    return tau_t
    
def dph_n(qh_t, ph_t, r, dr, d2r):
    #Lei de controle linear
    e = r - qh_t
    de = dr - ph_t
    #kp = 1600.0*eye(2)
    #kv = 2.0*sqrt(1600.0)*eye(2)
    #uh_lin_t = SMatrix(d2r + kv*de + kp*e, ph_.rowl_)
    lamda = 500.0
    k = 20.0
    n = 100.0
    Tanh = lambda x_: Matrix([tanh(x_[i]) for i in xrange(x_.rows)])
    uh_lin_t = SMatrix(d2r + lamda*de + k*Tanh(n*(de + lamda*e)), ph_.rowl_)        
    
    ph_t = SMatrix( ph_t, ph_.rowl_)
    q_t = fInvKin(qh_t)
    Ah_t = Ah_n(q_t)
    Ao_t = Ao_n(q_t)
    iAo_t = Ao_t.inv()
    C_t = SMatrix(1, ph_.rowl_, ph_.rowl_) - iAo_t*Ah_t
    p_t = C_t*ph_t
    b_t = b_n(q_t, p_t.M_)
    dCph_t = SMatrix(0, ph_.rowl_) + iAo_t*b_t
    M_t = M_n(q_t)
    v_t = v_n(q_t, p_t.M_)
    g_t = g_n(q_t)
    Mh_t = C_t.T()*M_t*C_t
    uh_t = Mh_t*uh_lin_t #+ C_t.T()*(M_t*dCph_t + v_t + g_t)
    u_t = SMatrix(0, p_.rowl_) + uh_t
    dph_t = Mh_t.LDLsolve( C_t.T()*( u_t - M_t*dCph_t - v_t - g_t ) )
    return dph_t
    
#def dph_n2(qh_t, ph_t, r, dr, d2r):
#    #Lei de controle linear
#    e = r - qh_t
#    de = dr - ph_t
#    #kp = 1600.0*eye(2)
#    #kv = 2.0*sqrt(1600.0)*eye(2)
#    #uh_lin_t = SMatrix(d2r + kv*de + kp*e, ph_.rowl_)
#    lamda = 500.0
#    k = 20.0
#    n = 100.0
#    Tanh = lambda x_: Matrix([tanh(x_[i]) for i in xrange(x_.rows)])
#    uh_lin_t = SMatrix(d2r + lamda*de + k*Tanh(n*(de + lamda*e)), ph_.rowl_)        
#    
#    ph_t = SMatrix( ph_t, ph_.rowl_)
#    q_t = fInvKin(qh_t)
#    Jh_t = Ah_n(q_t)
#    Jo_t = Ao_n(q_t)
#    C0_t = C0_n(q_t)
#    C1_t = C1_n(q_t)
#    C2_t = C2_n(q_t)
#    iJo_t = Jo_t.inv()
#    
#    Cch_t = SMatrix(1, ph_.rowl_, ph_.rowl_) - iJo_t*Jh_t
#    C_t = (C0_t + C1_t + C2_t)*Cch_t
#    p_t = C_t*ph_t
#    b1_t = b1_n(q_t, p_t.M_)
#    b2_t = b2_n(q_t, p_t.M_)
#    bch_t = bch_n(q_t, p_t.M_)
#    dCph_t = (C0_t + C1_t + C2_t)*( SMatrix(0, ph_.rowl_, ph_.coll_) + iJo_t*bch_t ) + (SMatrix(0, p_.rowl_, p_.coll_) - b1_t - b2_t )
#    M_t = M_n(q_t)
#    v_t = v_n(q_t, p_t.M_)
#    g_t = g_n(q_t)
#    Mh_t = C_t.T()*M_t*C_t
#    uh_t = Mh_t*uh_lin_t #+ C_t.T()*(M_t*dCph_t + v_t + g_t)
#    u_t = SMatrix(0, p_.rowl_) + uh_t
#    dph_t = Mh_t.LDLsolve( C_t.T()*( u_t - M_t*dCph_t - v_t - g_t ) )
#    return dph_t