import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from random import randint
from math import *
import numpy as np
from scipy.linalg import expm
from numpy.linalg import solve
from matplotlib import rc

rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

Eye = lambda N:np.matrix( np.eye(N) )
Zeros = lambda Ni,Nj :np.matrix( np.zeros((Ni,Nj)) )

class RP:
  def __init__(self, l1, l2, lg1, lg2, m1, m2, Jz1, Jz2, g):
    self.l1 = l1
    self.l2 = l2
    self.lg1 = lg1
    self.lg2 = lg2
    self.m1 = m1
    self.m2 = m2
    self.Jz1 = Jz1
    self.Jz2 = Jz2
    self.g = g
    
    self.q = Zeros(2,1)
    self.dq = Zeros(2,1)
    self.d2q = Zeros(2,1)
    
    self.x = Zeros(4,1)
    
    self.M = Zeros(2,2)
    self.V = Zeros(2,2)
    self.G = Zeros(2,2)
    self.h = Zeros(2,1)
    
    self.A11 = Zeros(2,2)
    self.A12 = Zeros(2,2)
    self.A21 = Eye(2)
    self.A22 = Zeros(2,2)
    
    self.B1 = Zeros(2,2)
    self.B2 = Zeros(2,2)
    
    self.A = Zeros(4,4)
    self.B = Zeros(4,2)
    self.f = Zeros(4,1)
    
  def doit(self, q, dq):
    self.q = q
    self.dq = dq
    
    self.x = np.vstack((self.dq, self.q))
    
    self.M = np.matrix([[self.Jz1+self.Jz2+self.m1*self.lg1**2 + self.m2*(q[1,0] -self.l2 + self.lg2 )**2, 0.0 ], [0.0,  self.m2]])
    
    self.h = np.matrix([2*self.m2*(q[1,0] -self.l2 + self.lg2)*dq[0,0]*dq[1,0]-(self.m1*self.lg1+self.m2*(dq[0,0] -self.l2 + self.lg2))*self.g*cos(q[0,0]), self.m2*(q[1,0] -self.l2 + self.lg2)*dq[0,0]**2 - self.m2*self.g*sin(q[0,0])]).T
    
  def doitx(self, x):
    self.q = x[2:4,0]
    self.dq = x[0:2,0]
    
    self.x = x
    
    q = self.q
    dq = self.dq
    d2q = self.d2q
    
    self.M = np.matrix([[self.Jz1+self.Jz2+self.m1*self.lg1**2 + self.m2*(q[1,0] -self.l2 + self.lg2 )**2, 0.0 ], [0.0,  self.m2]])
    
    self.h = np.matrix([2*self.m2*(q[1,0] -self.l2 + self.lg2)*dq[0,0]*dq[1,0]-(self.m1*self.lg1+self.m2*(dq[0,0] -self.l2 + self.lg2))*self.g*cos(q[0,0]), self.m2*(q[1,0] -self.l2 + self.lg2)*dq[0,0]**2 - self.m2*self.g*sin(q[0,0])]).T  
  
  def get_u(self, u):
    self.d2q = solve( self.M, (u - self.h) )
    self.f = np.vstack((self.d2q , self.dq  ))
    
    q = self.q
    dq = self.dq
    d2q = self.d2q
    
    self.V = np.matrix([[2*self.m2*(q[1,0] -self.l2 + self.lg2)*dq[1,0] , 2*self.m2*(q[1,0] -self.l2 + self.lg2)*dq[0,0] ], [-2*self.m2*(q[1,0] -self.l2 + self.lg2)*dq[0,0], 0.0]])
    
    self.G = np.matrix([[(self.m1*self.lg1+self.m2*(q[1,0] -self.l2 + self.lg2))*self.g*sin(q[0,0]), self.m2*(2*dq[1,0]+2*(q[1,0] -self.l2 + self.lg2)*d2q[0,0] - self.g*cos(q[0,0]) ) ], [-self.m2*self.g*cos(q[0,0]), -self.m2*d2q[0,0]]])
    
    self.A11 = solve( -(self.M), self.V )
    self.A12 = solve( -(self.M), self.G )
    self.B1 = (self.M).I
    
    self.A = np.vstack(( np.vstack((self.A11.T,self.A12.T)).T , np.vstack((self.A21.T,self.A22.T)).T ))
    
    self.B = np.vstack((self.B1, self.B2))
    
  
  
  def get_d2q(self, d2q):
    self.d2q = d2q
    self.f = np.vstack((self.d2q, self.dq))    
    
    q = self.q
    dq = self.dq
    
    
    self.V = np.matrix([[2*self.m2*(q[1,0] -self.l2 + self.lg2)*dq[1,0] , 2*self.m2*(q[1,0] -self.l2 + self.lg2)*dq[0,0] ], [-2*self.m2*(q[1,0] -self.l2 + self.lg2)*dq[0,0], 0.0]])
    
    self.G = np.matrix([[(self.m1*self.lg1+self.m2*(q[1,0] -self.l2 + self.lg2))*self.g*sin(q[0,0]), self.m2*(2*dq[1,0]+2*(q[1,0] -self.l2 + self.lg2)*d2q[0,0] - self.g*cos(q[0,0]) ) ], [-self.m2*self.g*cos(q[0,0]), -self.m2*d2q[0,0]]])
    
    self.A11 = solve(-(self.M), self.V)
    self.A12 = solve( -(self.M), self.G)
    self.B1 = (self.M).I
    
    self.A = np.vstack(( np.vstack((self.A11.T,self.A12.T)).T , np.vstack((self.A21.T,self.A22.T)).T ))
    
    self.B = np.vstack((self.B1, self.B2))
    
    
    
l1e = 0.1
l2e = 0.1
lg1e = 0.05
lg2e = 0.05
m1e = 0.1
m2e = 0.1
Jz1e = 80.0e-6
Jz2e = 80.0e-6
ge = 9.8

l1 = np.random.normal(l1e,0.5*0.001)
l2 = np.random.normal(l2e,0.5*0.001)
lg1 = np.random.normal(lg1e,0.5*0.002)
lg2 = np.random.normal(lg2e,0.5*0.002)
m1 = np.random.normal(m1e,0.5*0.1*m1e)
m2 = np.random.normal(m2e,0.5*0.1*m2e)
Jz1 = np.random.normal(Jz1e,0.5*0.2*Jz1e)
Jz2 = np.random.normal(Jz2e,0.5*0.2*Jz2e)
g = ge
    
R1 = RP(l1=l1, l2=l2, lg1=lg1, lg2=lg2, m1=m1, m2=m2, Jz1=Jz1, Jz2=Jz2, g=g)
R2 = RP(l1=l1, l2=l2, lg1=lg1, lg2=lg2, m1=m1, m2=m2, Jz1=Jz1, Jz2=Jz2, g=g)

w = 10.0
theta_d = lambda t: 1.0*(-pi/2) + 0.5*pi*sin(w*t)
dtheta_d = lambda t:w*0.5*pi*cos(w*t)
d2theta_d = lambda t:-w*w*0.5*pi*sin(w*t)
d_d = lambda t:0.15 + 0.1*sin(w*t)
dd_d = lambda t:w*0.1*cos(w*t)
d2d_d = lambda t:-w*w*0.1*sin(w*t)

qd = lambda t:np.matrix([theta_d(t),d_d(t)]).T
dqd = lambda t:np.matrix([dtheta_d(t),dd_d(t)]).T
d2qd = lambda t:np.matrix([d2theta_d(t),d2d_d(t)]).T

T = 5e-6
tf = 0.2*6.28
n = int(round(tf/T + 1) )
t_ = [k*T for k in range(n)]
u_ = [Zeros(2,1) for t in t_]
ud_ = [Zeros(2,1) for t in t_]
x_ = [Zeros(4,1) for t in t_]
z_ = [Zeros(4,1) for t in t_]
xe_ = [Zeros(4,1) for t in t_]
P_ = [Zeros(4,4) for t in t_]

mean = 0
std1 = 0.5*2.0*pi/4000.0
std2 = 0.5*0.001
num_samples = n
v1 = np.random.normal(mean, std1, size=num_samples)
v2 = np.random.normal(mean, std2, size=num_samples)
v_= [np.matrix([v1[k],v2[k]]).T for k in range(n)]

R = np.matrix([[std1**2, 0],[0,std2**2]])
Q = np.matrix([[1*std1**2,0,0,0],[0,1*std2**2,0,0],[0,0,0,0],[0,0,0,0]])

H = np.vstack((Zeros(2,2),Eye(2))).T  

Ad = Zeros(4,4)
K1 = Zeros(4,1)
K2 = Zeros(4,1)
K3 = Zeros(4,1)
K4 = Zeros(4,1)

lamb = 100.0
kp = lamb**2
kv = 2*lamb

R1.doit(qd(0),Zeros(2,1))
R2.doit(qd(0),dqd(0))
R2.get_d2q(d2qd(0))
x_[0] = R1.x
xe_[0] = R2.x 
P_[0] = np.matrix([[(0.5*pi*w)**2, 0,0,0],[0,(0.1*w)**2,0,0],[0,0,std1**2,0],[0,0,0,std2**2]])
  
z_[0] = H*x_[0] + v_[0]
#K =  P_[0]*H.T*( H*P_[0]*H.T + R ).I
K = solve( (H*P_[0]*H.T + R).T , H*P_[0].T ).T
xe_[0] = xe_[0] + K*(z_[0] - H*xe_[0])
P_[0] = (Eye(4) - K*H)*P_[0]*(Eye(4) - K*H).T + K*R*K.T

R2.doitx(xe_[0])

for k in range(n-1):
  u_[k] = R2.h +  R2.M*(d2qd(t_[k]) + kv*(dqd(t_[k]) - R2.dq ) + kp*(qd(t_[k]) - R2.q )  )
  R1.get_u(u_[k])
  R2.get_u(u_[k])

  K1 = R1.f
  R1.doitx(x_[k]+ 0.5*T*K1)
  R1.get_u(u_[k])
  K2 = R1.f
  R1.doitx(x_[k]+ 0.5*T*K2)
  R1.get_u(u_[k])  
  K3 = R1.f
  R1.doitx(x_[k]+ T*K3)
  R1.get_u(u_[k])
  K4 = R1.f
  x_[k+1] = x_[k] + (T/6.0)*(K1 + 2*K2 + 2*K3 + K4)
  R1.doitx(x_[k+1])
  """
  R1.doitx(x_[k]+ 0.5*T*K1)
  R1.get_u(u_[k])
  K2 = R1.f
  R1.doitx(x_[k]+ T*((2.0/9)*K1 + (4.0/9)*K2 ) )
  R1.get_u(u_[k])  
  K3 = R1.f
  R1.doitx(x_[k]+ T*((7.0/36)*K1 + (2.0/9)*K2 + (-1.0/12)*K3 )  )
  R1.get_u(u_[k])
  K4 = R1.f
  R1.doitx(x_[k]+ T*((-35.0/144)*K1 + (-55.0/36)*K2 + (35.0/48)*K3 + (15.0/8)*K4 )  )
  R1.get_u(u_[k])
  K5 = R1.f
  R1.doitx(x_[k]+ T*((-1.0/360)*K1 + (-11.0/36)*K2 + (-1.0/8)*K3 + (0.5)*K4 + (1.0/10)*K5 )  )
  R1.get_u(u_[k])
  K6 = R1.f
  R1.doitx(x_[k]+ T*((-41.0/260)*K1 + (22.0/13)*K2 + (43.0/156)*K3 + (-118.0/39)*K4 + (32.0/195)*K5 + (80.0/39)*K6 )  )
  R1.get_u(u_[k])
  K7 = R1.f
  x_[k+1] = x_[k] + T*( (13.0/200)*K1 + (11.0/40)*K3 + (11.0/40)*K4 + (4.0/25)*K5 + (4.0/25)*K6 + (13.0/200)*K7 ) 
  R1.doitx(x_[k+1])
  """
  
  #Predicao
  Ad = Eye(4) + T*R2.A + (T*R2.A**2)**2/2 + (T*R2.A)**3/6 +  (T*R2.A)**4/24
  #Ad = Eye(4) + T*R2.A + (T*R2.A**2)**2/2 + (T*R2.A)**3/6 +  (T*R2.A)**4/24 + (T*R2.A)**5/120 + (T*R2.A)**6/720
  #Ad = np.matrix( expm(R2.A) )
  P_[k+1] = Ad*P_[k]*Ad.T + Q
  
  K1 = R2.f
  R2.doitx(xe_[k]+ 0.5*T*K1)
  R2.get_u(u_[k])
  K2 = R2.f
  R2.doitx(xe_[k]+ 0.5*T*K2)
  R2.get_u(u_[k])  
  K3 = R2.f
  R2.doitx(xe_[k]+ T*K3)
  R2.get_u(u_[k])
  K4 = R2.f
  xe_[k+1] = xe_[k] + (T/6.0)*(K1 + 2*K2 + 2*K3 + K4)
  """
  R2.doitx(xe_[k]+ 0.5*T*K1)
  R2.get_u(u_[k])
  K2 = R2.f
  R2.doitx(xe_[k]+ T*((2.0/9)*K1 + (4.0/9)*K2 ) )
  R2.get_u(u_[k])  
  K3 = R2.f
  R2.doitx(xe_[k]+ T*((7.0/36)*K1 + (2.0/9)*K2 + (-1.0/12)*K3 )  )
  R2.get_u(u_[k])
  K4 = R2.f
  R2.doitx(xe_[k]+ T*((-35.0/144)*K1 + (-55.0/36)*K2 + (35.0/48)*K3 + (15.0/8)*K4 )  )
  R2.get_u(u_[k])
  K5 = R2.f
  R2.doitx(xe_[k]+ T*((-1.0/360)*K1 + (-11.0/36)*K2 + (-1.0/8)*K3 + (0.5)*K4 + (1.0/10)*K5 )  )
  R2.get_u(u_[k])
  K6 = R2.f
  R2.doitx(xe_[k]+ T*((-41.0/260)*K1 + (22.0/13)*K2 + (43.0/156)*K3 + (-118.0/39)*K4 + (32.0/195)*K5 + (80.0/39)*K6 )  )
  R2.get_u(u_[k])
  K7 = R2.f
  xe_[k+1] = xe_[k] + T*( (13.0/200)*K1 + (11.0/40)*K3 + (11.0/40)*K4 + (4.0/25)*K5 + (4.0/25)*K6 + (13.0/200)*K7 )
  """
  #Correcao
  z_[k+1] = H*x_[k+1] + v_[k+1]
  #K =  P_[k+1]*H.T * ( H*P_[k+1]*H.T + R ).I
  K = solve( (H*P_[k+1]*H.T + R).T , H*P_[k+1].T ).T
  xe_[k+1] = xe_[k+1] + K*(z_[k+1] - H*xe_[k+1])
  P_[k+1] = (Eye(4) - K*H)*P_[k+1]*(Eye(4) - K*H).T + K*R*K.T
  #Atualiza modelo
  R2.doitx(xe_[k+1])
  print(k)
  

t_np = np.array([t_[i] for i in range(n)])
theta_np = np.array([x_[i][2,0] for i in range(n)])
d_np = np.array([x_[i][3,0] for i in range(n)])
dtheta_np = np.array([x_[i][0,0] for i in range(n)])
dd_np = np.array([x_[i][1,0] for i in range(n)])
thetae_np = np.array([xe_[i][2,0] for i in range(n)])
de_np = np.array([xe_[i][3,0] for i in range(n)])
dthetae_np = np.array([xe_[i][0,0] for i in range(n)])
dde_np = np.array([xe_[i][1,0] for i in range(n)])
theta_d_np = np.array([qd(t_[i])[0,0] for i in range(n)])
d_d_np = np.array([qd(t_[i])[1,0] for i in range(n)])
dtheta_d_np = np.array([dqd(t_[i])[0,0] for i in range(n)])
dd_d_np = np.array([dqd(t_[i])[1,0] for i in range(n)])
e_theta_np = theta_d_np - theta_np
e_d_np = d_d_np - d_np

e_obs_theta_np = theta_np - thetae_np
e_obs_d_np = d_np - de_np
e_obs_dtheta_np = dtheta_np - dthetae_np
e_obs_dd_np = dd_np - dde_np

P00_np = np.array([P_[i][0,0] for i in range(n)])
P11_np = np.array([P_[i][1,1] for i in range(n)])
P22_np = np.array([P_[i][2,2] for i in range(n)])
P33_np = np.array([P_[i][3,3] for i in range(n)])

Ptr_np = P00_np + P11_np + P22_np + P33_np


fig1 = plt.figure()
plt.plot(t_np, theta_np, 'r', t_np, thetae_np, 'k--', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$\theta [rad]$')
#plt.title('s em z')
plt.savefig('theta.pdf')

plt.figure()
plt.plot(t_np, d_np, 'r', t_np, de_np, 'k--', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$d[m]$')
#plt.title('s em z')
plt.savefig('d.pdf')

plt.figure()
plt.plot(t_np, dtheta_np, 'r', t_np, dthetae_np, 'k--', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$\dot{\theta}[rad/s]$')
#plt.title('s em z')
plt.savefig('dtheta.pdf')

plt.figure()
plt.plot(t_np, dd_np, 'r', t_np, dde_np, 'k--', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$\dot{d}[m/s]$')
#plt.title('s em z')
plt.savefig('dd.pdf')

plt.figure()
plt.plot(t_np, e_theta_np, 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$ (\theta_d - \theta) [rad]  $')
#plt.title('s em z')
plt.savefig('e_theta.pdf')

plt.figure()
plt.plot(t_np, e_d_np, 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$ (d_d - d) [m]  $')
#plt.title('s em z')
plt.savefig('e_d.pdf')

plt.figure()
plt.plot(t_np[0:int(n/1000.0)], P00_np[0:int(n/1000.0)], 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$P_{\dot{\theta},\dot{\theta}}$')
#plt.title('s em z')
plt.savefig('P00.pdf')

plt.figure()
plt.plot(t_np[0:int(n/1000.0)], P11_np[0:int(n/1000.0)], 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$P_{\dot{d},\dot{d}}$')
#plt.title('s em z')
plt.savefig('P11.pdf')

plt.figure()
plt.plot(t_np[0:int(n/1000.0)], P22_np[0:int(n/1000.0)], 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$P_{\theta,\theta}$')
#plt.title('s em z')
plt.savefig('P22.pdf')

plt.figure()
plt.plot(t_np[0:int(n/1000.0)], P33_np[0:int(n/1000.0)], 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$P_{d,d}$')
#plt.title('s em z')
plt.savefig('P33.pdf')

plt.figure()
plt.plot(t_np[0:int(n/1000.0)], Ptr_np[0:int(n/1000.0)], 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$tr(\mathbf{P})$')
#plt.title('s em z')
plt.savefig('Ptr.pdf')

plt.figure()
plt.plot(t_np[0:int(n/8.3)], e_obs_theta_np[0:int(n/8.3)], 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$( \theta - \hat{\theta} )[rad]$')
#plt.title('s em z')
plt.savefig('e_obs_theta.pdf')

plt.figure()
plt.plot(t_np[0:int(n/8.3)], e_obs_d_np[0:int(n/8.3)], 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$( \theta - \hat{\theta} )[rad]$')
#plt.title('s em z')
plt.savefig('e_obs_d.pdf')

plt.figure()
plt.plot(t_np[0:int(n/500.0)], e_obs_dtheta_np[0:int(n/500.0)], 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$( \dot{\theta} - \dot{\hat{\theta}} )[rad/s]$')
#plt.title('s em z')
plt.savefig('e_obs_dtheta.pdf')

plt.figure()
plt.plot(t_np[0:int(n/500.0)], e_obs_dd_np[0:int(n/500.0)], 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$( \dot{d} - \dot{\hat{d}} )[m/s]$')
#plt.title('s em z')
plt.savefig('e_obs_dd.pdf')

