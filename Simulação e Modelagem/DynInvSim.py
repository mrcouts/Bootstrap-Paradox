from DynInvNum import *
import matplotlib.pyplot as plt
import numpy as np
from math import pi as Pi
   
#Parametros da trajetoria
r = 0.085
x0 = 0.0
y0 = 0.16

v=1.0 #Velocidade tangencial
w=v/r #velocidade angular
p=0.6 #tempo total da simulacao
dt=0.003 #dt
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