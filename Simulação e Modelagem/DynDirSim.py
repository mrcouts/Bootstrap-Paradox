from DynInvNum import *
from RK import *
from math import pi as Pi
   
#Parametros da trajetoria
r = 0.04
x0 = 0.0
y0 = 0.16

v=0.4 #Velocidade tangencial
w=v/r #velocidade angular
p=0.6 #tempo total da simulacao
dt=0.001 #dt
#nt=int(p/dt)+1 #tamanho do vetor de tempo
#t_=[i*dt for i in xrange(nt)] #vetor de tempo

#Trajetoria
r_ = lambda T: Matrix([x0+r*cos(w*T), y0+r*sin(w*T)])

#velocidades
dr_ = lambda T: Matrix([-w*r*sin(w*T), w*r*cos(w*T)])

#aceleracoes
d2r_ = lambda T: Matrix([-(w**2)*r*cos(w*T), -(w**2)*r*sin(w*T)])

#Runge-Kutta
y0_ = Matrix([r_(0.0),dr_(0.0)])
f = lambda T, y_: Matrix([ y_[2:4, 0], dph_n( y_[0:2, 0] , y_[2:4, 0], r_(T), dr_(T), d2r_(T)).M_ ])

h=dt
tf=p
rk=RK()
y_,t_=rk.aplic(h,tf,y0_,f)

#Plotando os graficos dos erros

import matplotlib.pyplot as plt
import numpy as np

t_np = np.linspace(t[0], t[-1], len(t))
ex_np = t_np.copy()
ey_np = t_np.copy()

for i in np.arange(np.size(t_np)):
    ex_np[i] = r_(t_[i])[0] - y_[i][0]
    ey_np[i] = r_(t_[i])[1] - y_[i][1]
    

plt.figure()
plt.plot(t_np, ex_np, 'r')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()

plt.figure()
plt.plot(t_np, ey_np, 'r')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()


"""tau_t = [tau_n(Xt[i], dXt[i], d2Xt[i]) for i in range(nt)]

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
"""