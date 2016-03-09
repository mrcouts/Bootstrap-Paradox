from math import pi as Pi
import matplotlib.pyplot as plt
import numpy as np
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

C_s = C_.subs(PARAMETROS)
Z_s = Z_.subs(PARAMETROS)
dC_s = dC_.subs(PARAMETROS)
M_s = M_.subs(PARAMETROS)
v_s = v_.subs(PARAMETROS)
g_s = g_.subs(PARAMETROS)

C_n = lambda q_n:       SMatrix( ( C_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), C_.rowl_, C_.coll_)
Z_n = lambda q_n:       SMatrix( ( Z_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), Z_.rowl_, Z_.coll_) 
M_n = lambda q_n:       SMatrix( ( M_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), M_.rowl_, M_.coll_)
g_n = lambda q_n:       SMatrix( ( g_s.M_.subs(rep( list(q_.M_) ,q_n))).evalf(), g_.rowl_, g_.coll_)
dC_n= lambda q_n, p_n:  SMatrix( (dC_s.M_.subs(rep( list(p_.M_), p_n)).subs(rep( list(q_.M_) ,q_n))).evalf(), dC_.rowl_, dC_.coll_)
v_n = lambda q_n, p_n:  SMatrix( ( v_s.M_.subs(rep( list(p_.M_), p_n)).subs(rep( list(q_.M_) ,q_n))).evalf(), v_.rowl_, v_.coll_)
  
def tau_n(qh_t, ph_t, dph_t):
    q_t = fInvKin(qh_t)
    C_t = C_n(q_t)
    p_t = C_t*SMatrix( ph_t, ph_.rowl_)
    dC_t = dC_n(q_t, p_t.M_)
    dp_t = dC_t*SMatrix( ph_t ,ph_.rowl_) + C_t*SMatrix( dph_t,ph_.rowl_)
    Z_t = Z_n(q_t)
    M_t = M_n(q_t)
    v_t = v_n(q_t, p_t.M_)
    g_t = g_n(q_t)    
    tau_t = Z_t.inv()*C_t.T()*( M_t*dp_t + v_t + g_t )
    return tau_t
    
#Parametros da trajetoria
r = 0.085
x0 = 0.0
y0 = 0.16

v=1.0 #Velocidade tangencial
w=v/r #velocidade angular
p=3*0.6 #tempo total da simulacao
dt=0.01 #dt
nt=int(p/dt)+1 #tamanho do vetor de tempo
t=[i*dt for i in xrange(nt)] #vetor de tempo

#Trajetoria
xt=[x0+r*cos(w*T) for T in t]
yt=[y0+r*sin(w*T) for T in t]
Xt=[Matrix([xt[i],yt[i]]) for i in xrange(nt)]

#velocidades
dxt=[-w*r*sin(w*T) for T in t]
dyt=[w*r*cos(w*T) for T in t]
dXt=[Matrix([dxt[i],dyt[i]]) for i in xrange(nt)]

#aceleracoes
d2xt=[-(w**2)*r*cos(w*T) for T in t]
d2yt=[-(w**2)*r*sin(w*T) for T in t]
d2Xt=[Matrix([d2xt[i],d2yt[i]]) for i in xrange(nt)]

tau_t = [tau_n(Xt[i], dXt[i], d2Xt[i]) for i in range(nt)]

t_np = np.linspace(t[0], t[-1], len(t))
tau1_np = t_np.copy()
tau2_np = t_np.copy()

for i in np.arange(np.size(t_np)):
    tau1_np[i] = tau_t[i].M_[0]
    tau2_np[i] = tau_t[i].M_[1]

plt.figure()
plt.plot(t_np, tau1_np, 'r')
#plt.xscale('log')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()

plt.figure()
plt.plot(t_np, tau2_np, 'r')
#plt.xscale('log')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()