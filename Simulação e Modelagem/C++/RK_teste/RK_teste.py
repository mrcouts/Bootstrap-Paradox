import matplotlib.pyplot as plt
import numpy as np
from math import pi as Pi
from txt2py import *

A = txt2py("RK_teste.txt")

t_np = np.array([A[i][0] for i in range(len(A))])
tau1_np = np.array([A[i][1] for i in range(len(A))])
tau2_np = np.array([A[i][2] for i in range(len(A))])

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