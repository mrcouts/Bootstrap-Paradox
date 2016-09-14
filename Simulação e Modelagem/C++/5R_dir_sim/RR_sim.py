import matplotlib.pyplot as plt
import numpy as np
from math import pi as Pi
from txt2py import *

A = txt2py("RR_sim.txt")

t_np = np.array([A[i][0] for i in range(len(A))])
tau1_np = np.array([A[i][1] for i in range(len(A))])
tau2_np = np.array([A[i][2] for i in range(len(A))])
erro1_np= np.array([A[i][3] for i in range(len(A))])
erro2_np= np.array([A[i][4] for i in range(len(A))])
s1_np   = np.array([A[i][5] for i in range(len(A))])
s2_np   = np.array([A[i][6] for i in range(len(A))])

plt.figure()
plt.plot(t_np, tau1_np, 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel('t[s]')
plt.ylabel('tau 1 [Nm]')
plt.title('Torque no atuador 1')
plt.show()

plt.figure()
plt.plot(t_np, tau2_np, 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel('t[s]')
plt.ylabel('tau 2 [Nm]')
plt.title('Torque no atuador 2')
plt.show()

plt.figure()
plt.plot(t_np, erro1_np, 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel('t[s]')
plt.ylabel('e1[m]')
plt.title('Erro em x')
plt.show()

plt.figure()
plt.plot(t_np, erro2_np, 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel('t[s]')
plt.ylabel('e2[m]')
plt.title('Erro em y')
plt.show()

plt.figure()
plt.plot(t_np, s1_np, 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel('t[s]')
plt.ylabel('s1[m/s]')
plt.title('s em x')
plt.show()

plt.figure()
plt.plot(t_np, s2_np, 'r', linewidth=2)
#plt.xscale('log')
plt.xlabel('t[s]')
plt.ylabel('s2[m/s]')
plt.title('s em y')
plt.show()