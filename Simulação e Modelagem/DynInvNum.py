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

Ah_n= lambda q_n:       SMatrix( (Ah_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), Ah_.rowl_, Ah_.coll_)
Ao_n= lambda q_n:       SMatrix( (Ao_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), Ao_.rowl_, Ao_.coll_)
M_n = lambda q_n:       SMatrix( ( M_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), M_.rowl_, M_.coll_)
g_n = lambda q_n:       SMatrix( ( g_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), g_.rowl_, g_.coll_)
v_n = lambda q_n, p_n:  SMatrix( ( v_s.M_.subs(rep( list(p_.M_), p_n)).subs(rep( list(q_.M_) ,q_n))).evalf(),  v_.rowl_,  v_.coll_)
b_n = lambda q_n, p_n:  SMatrix( ( b_s.M_.subs(rep( list(p_.M_), p_n)).subs(rep( list(q_.M_) ,q_n))).evalf(),  b_.rowl_,  b_.coll_)
      
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
    kp = 800.0*eye(2)
    kv = 2.0*sqrt(800.0)*eye(2)
    uh_lin_t = SMatrix(d2r + kv*de + kp*e, uh_.rowl_)        
    
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
    uh_t = Mh_t*uh_lin_t + C_t.T()*(M_t*dCph_t + v_t + g_t)
    u_t = SMatrix(0, p_.rowl_) + uh_t
    dph_t = Mh_t.LDLsolve( C_t.T()*( u_t - M_t*dCph_t - v_t - g_t ) )
    return dph_t