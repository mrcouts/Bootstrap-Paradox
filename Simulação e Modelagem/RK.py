from sympy import *
from math import pi as Pi
from numpy.linalg.linalg import dot
init_printing(use_unicode=True)

class RK (object):
    def __init__(self, method='RK4'):
        if method=='RK4':
            self.a=[[],[0.5],[0,0.5],[0,0,1.0]]
            self.b=[1.0/6, 1.0/3, 1.0/3, 1.0/6]
            self.c=[0.0]+[sum(x) for x in self.a]
            self.n=len(self.b)
            
        else:
            raise ValueError('Este metodo nao esta disponivel')
    
    def aplic(self,h,tf,y0,f):
        nt=int(tf/h)
        t=[x*h for x in xrange(nt+1)]
        y=[y0] + [0 for i in xrange(nt)]
        k=[0 for i in xrange(self.n)]         
        
        for i in xrange(nt):
            k[0]=f(t[i],y[i])
            for j in xrange(1,self.n):
                k[j]=f(t[i]+self.c[j]*h, y[i]+h*dot(self.a[j], k[0:j]))
            y[i+1]=y[i]+h*dot(self.b,k)
        
        return y,t

## f=lambda t,y: 1000.0*(sign(sin(1000.0*t))-y)
## y0=0.0
wn=1000.0
z=0.5
f=lambda t,y: Matrix([y[1], (wn**2)*(sign(sin(1000.0*t))-y[0])-2*wn*z*y[1]])
y0 = zeros(2,1)
h=1.0e-4
tf=4.0e-2
v=RK()
y,t=v.aplic(h,tf,y0,f)

import matplotlib.pyplot as plt
import numpy as np

t_np = np.linspace(t[0], t[-1], len(t))
y_np = t_np.copy()

for i in np.arange(np.size(t_np)):
    y_np[i] = y[i][0]
    

plt.figure()
plt.plot(t_np, y_np, 'r')
#plt.xscale('log')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()
