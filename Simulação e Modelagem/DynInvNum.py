from math import pi as Pi
import matplotlib.pyplot as plt
import numpy as np
from InKinNum import *
from DynamicsDeduction import *

PARAMETROS = [(symbols('l_0'), par[l0]),
              (symbols('l_1'), par[l1]),
              (symbols('l_2'), par[l2]),
              (symbols('lg_1'), 0.5*par[l1]),
              (symbols('lg_2'), 0.5*par[l2]),
              (symbols('m_1'), 0.143),
              (symbols('m_2'), 0.171),
              (symbols('Jz_1'), (0.143*par[l1]**2)/12.0),
              (symbols('Jz_2'), (0.171*par[l2]**2)/12.0),
              (symbols('g'), 9.8)         
]
                  
rep = lambda vec1,vec2: [(vec1[i],vec2[i]) for i in xrange(len(vec1))]
    
C_n = lambda q_n: SMatrix( (C_.M_.subs(rep( list(q_.M_) ,q_n)).subs(PARAMETROS)).evalf(), C_.rowl_, C_.coll_)
dC_n = lambda q_n, p_n: SMatrix( (dC_.M_.subs( rep(list(p_.M_), p_n) ).subs(rep( list(q_.M_) ,q_n)).subs(PARAMETROS)).evalf(), dC_.rowl_, dC_.coll_)  
Z_n = lambda q_n: SMatrix( (Z_.M_.subs(rep( list(q_.M_) ,q_n)).subs(PARAMETROS)).evalf(), Z_.rowl_, Z_.coll_)  


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

q_t = []
C_t = []
dC_t = []
Z_t = []
p_t = []
dp_t = []
tau_t = []

M_n = SMatrix( M_.M_.subs(PARAMETROS), M_.rowl_, M_.coll_)
v_n = SMatrix( v_.M_.subs(PARAMETROS), v_.rowl_)
g_n = SMatrix( g_.M_.subs(PARAMETROS), g_.rowl_)

for i in xrange(nt):
    q_t.append(fInvKin(Xt[i]))
    C_t.append(C_n(q_t[i]))
    Z_t.append(Z_n(q_t[i]))
    p_t.append(C_t[i]*SMatrix( dXt[i],ph_.rowl_))
    dC_t.append(dC_n(q_t[i], p_t[i].M_))
    dp_t.append( dC_t[i]*SMatrix( dXt[i],ph_.rowl_) + C_t[i]*SMatrix( d2Xt[i],ph_.rowl_) )
    tau_t.append( Z_t[i].inv()*C_t[i].T()*( M_n*dp_t[i] + v_n + g_n ) )

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