from DynInvNum import *
from RK import *
from math import pi as Pi
   
#Parametros da trajetoria
r = 0.08
x0 = 0.0
y0 = 0.16

v=1.0 #Velocidade tangencial
w=v/r #velocidade angular

#Trajetoria
r_ = lambda T: Matrix([x0+r*cos(w*T), y0+r*sin(w*T)])

#velocidades
dr_ = lambda T: Matrix([-w*r*sin(w*T), w*r*cos(w*T)])

#aceleracoes
d2r_ = lambda T: Matrix([-(w**2)*r*cos(w*T), -(w**2)*r*sin(w*T)])

#Runge-Kutta
y0_ = Matrix([r_(0.0),0.0*dr_(0.0)])

def f (T, y_):
     dph_t, uh_t = dph_n( y_[0:2, 0] , y_[2:4, 0], r_(T), dr_(T), d2r_(T))
     return Matrix([ y_[2:4, 0], dph_t.M_ ]), uh_t.M_

h=0.006 #dt
tf=2*0.6 #tempo total da simulacao
rk=RK('RK6')

import time
start = time.time()
y_,u_t,t_=rk.Apply2(h,tf,y0_,f)
end = time.time()
print(end - start)

#Plotando os graficos dos erros

import matplotlib.pyplot as plt
import numpy as np

t_np = np.linspace(t_[0], t_[-1], len(t_))
ex_np = t_np.copy()
ey_np = t_np.copy()
u1_np = t_np.copy()
u2_np = t_np.copy()

for i in np.arange(np.size(t_np)):
    ex_np[i] = r_(t_[i])[0] - y_[i][0]
    ey_np[i] = r_(t_[i])[1] - y_[i][1]
    u1_np[i] = u_t[i][0]
    u2_np[i] = u_t[i][1]
    

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

plt.figure()
plt.plot(t_np, u1_np, 'r')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()

plt.figure()
plt.plot(t_np, u2_np, 'r')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()