# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from math import pi as Pi
from txt2py import *

A = txt2py("RR_sim.txt")

t_np = np.array([A[i][0] for i in range(len(A))])
tau1_np = np.array([A[i][1] for i in range(len(A))])
tau2_np = np.array([A[i][2] for i in range(len(A))])
erro1_np= 1000.0*np.array([A[i][3] for i in range(len(A))])
erro2_np= 1000.0*np.array([A[i][4] for i in range(len(A))])
s1_np   = np.array([A[i][5] for i in range(len(A))])
s2_np   = np.array([A[i][6] for i in range(len(A))])

plt.figure()
plt.plot(t_np, tau1_np, 'r', linewidth=2)
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$\tau_1[Nm]$')
plt.ylim(-0.2, 1.0) 
plt.title('Torque aplicado ' + 'em fun'  + u'ç' + u'ã' 'o do tempo (atuador 1)')
plt.savefig('tau1.png')

plt.figure()
plt.plot(t_np, tau2_np, 'r', linewidth=2)
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$\tau_2[Nm]$')
plt.ylim(-0.2, 0.4)
plt.title('Torque aplicado ' + 'em fun'  + u'ç' + u'ã' 'o do tempo (atuador 2)')
plt.savefig('tau2.png')

plt.figure()
plt.plot(t_np, erro1_np, 'r', linewidth=2)
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$e_x[mm]$')
plt.ylim(-4, 3)
plt.title('Erro de posi' + u'ç' + u'ã' + 'o em fun'  + u'ç' + u'ã' 'o do tempo (coordenada x)')
plt.savefig('ex.png')

plt.figure()
plt.plot(t_np, erro2_np, 'r', linewidth=2)
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$e_y[mm]$')
plt.ylim(-6, 6)
plt.title('Erro de posi' + u'ç' + u'ã' + 'o em fun'  + u'ç' + u'ã' 'o do tempo (coordenada y)')
plt.savefig('ey.png')

#plt.figure()
#plt.plot(t_np, s1_np, 'r', linewidth=2)
##plt.xscale('log')
#plt.xlabel('t[s]')
#plt.ylabel('s1[m/s]')
#plt.title('s em x')
#plt.show()
#
#plt.figure()
#plt.plot(t_np, s2_np, 'r', linewidth=2)
##plt.xscale('log')
#plt.xlabel('t[s]')
#plt.ylabel('s2[m/s]')
#plt.title('s em y')
#plt.show()