from sympy import *
from math import pi as Pi
from numpy.linalg.linalg import dot
init_printing(use_unicode=True)

class RK(object):
    def __init__(self, method='RK6'):
        if method == 'Heun':
            self.a = [[1.0]]
            self.b = [0.5, 0.5]
        elif method == 'RK3':
            self.a = [[0.5],
                      [-1.0, 2.0]]
            self.b = [1.0/6, 2.0/3, 1.0/6]
        elif method=='RK4':
            self.a=[[0.5],
                    [0,0.5],
                    [0,0,1.0]]
            self.b=[1.0/6, 1.0/3, 1.0/3, 1.0/6]
        elif method=='RK6':
            self.a=[[0.5],
                    [2.0/9, 4.0/9],
                    [7.0/36, 2.0/9, -1.0/12],
                    [-35.0/144, -55.0/36, 35.0/48, 15.0/8],
                    [-1.0/360, -11.0/36, -1.0/8, 0.5, 1.0/10],
                    [-41.0/260, 22.0/13, 43.0/156, -118.0/39, 32.0/195, 80.0/39]]
            self.b=[13.0/200, 0, 11.0/40, 11.0/40, 4.0/25, 4.0/25, 13.0/200]
        else:
            raise ValueError('Este metodo nao esta disponivel')
        self.c=[sum(x) for x in self.a]
        self.n=len(self.b)
     
    def Apply(self,h,tf,y0_,f_):
        nt=int(tf/h)
        t_ = [i*h for i in xrange(nt+1)]
        y__= [y0_ for i in xrange(nt+1)]
        k__= [y0_ for i in xrange(self.n)]         
        
        for i in xrange(nt):
            k__[0] = f_(t_[i],y__[i])
            for j in xrange(1, self.n):
                k__[j] = f_(t_[i] + self.c[j-1]*h, y__[i] + h*dot(self.a[j-1], k__[0:j]))
            y__[i+1] = y__[i] + h*dot(self.b, k__)
            print(i)
        return y__, t_
                
    def Apply2(self,h,tf,y0_,f_):
        nt=int(tf/h)
        t_ = [i*h for i in xrange(nt+1)]
        y__= [y0_ for i in xrange(nt+1)]
        k__= [y0_ for i in xrange(self.n)]
        u__= [None for i in xrange(nt+1)]         
        
        for i in xrange(nt):
            k__[0], u__[i] = f_(t_[i], y__[i])
            for j in xrange(1, self.n):
                k__[j] = f_(t_[i] + self.c[j-1]*h, y__[i] + h*dot(self.a[j-1], k__[0:j]))[0]
            y__[i+1] = y__[i] + h*dot(self.b, k__)
            print(i)
        u__[nt] = f_(t_[nt], y__[nt])[1]
        return y__, u__, t_
        
## f=lambda t,y: 1000.0*(sign(sin(1000.0*t))-y)
## y0=0.0

"""wn=1000.0
z=0.5
f=lambda t,y: Matrix([y[1], (wn**2)*(sign(sin(1000.0*t))-y[0])-2*wn*z*y[1]])
y0 = zeros(2,1)
h=1.0e-4
tf=4.0e-2
v=RK()
y,t=v.Apply(h,tf,y0,f)

import matplotlib.pyplot as plt
import numpy as np

t_np = np.linspace(t[0], t[-1], len(t))
y_np = t_np.copy()

for i in np.arange(np.size(t_np)):
    y_np[i] = y[i][0]
    

plt.figure()
plt.plot(t_np, y_np, 'r')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()"""
