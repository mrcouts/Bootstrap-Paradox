import matplotlib.pyplot as plt
import numpy as np
from math import pi as Pi

def txt2py(filename):
    with open(filename, 'r') as MyFile:
        A = []
        cont = 0
        for line in MyFile:
            A.append(line.split(';'))
            A[cont].pop()
            cont += 1
    return [[float(A[i][j]) for j in range(len(A[0]))] for i in range(len(A))]

A = txt2py("6Rsim.txt")
print np.array([A[i][0] for i in range(len(A))])

t_np = np.array([A[i][0] for i in range(len(A))])
tau1_np = np.array([A[i][1] for i in range(len(A))])
tau2_np = np.array([A[i][2] for i in range(len(A))])
tau3_np = np.array([A[i][3] for i in range(len(A))])
tau4_np = np.array([A[i][4] for i in range(len(A))])
tau5_np = np.array([A[i][5] for i in range(len(A))])
tau6_np = np.array([A[i][6] for i in range(len(A))])

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

plt.figure()
plt.plot(t_np, tau3_np, 'r')
#plt.xscale('log')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()

plt.figure()
plt.plot(t_np, tau4_np, 'r')
#plt.xscale('log')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()

plt.figure()
plt.plot(t_np, tau5_np, 'r')
#plt.xscale('log')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()

plt.figure()
plt.plot(t_np, tau6_np, 'r')
#plt.xscale('log')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()