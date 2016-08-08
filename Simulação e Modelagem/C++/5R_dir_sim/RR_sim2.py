import matplotlib.pyplot as plt
import numpy as np
from math import pi as Pi
from txt2py import *
from sympy import *

A = txt2py("RR_sim2.txt")

rows = len(A)
cols = len(A[0])

M = zeros(rows, cols)
for i in range(rows):
    for j in range(cols):
        M[i,j] = A[i][j]

M2=M.extract(range(rows-1,-1,-1),range(cols))

M3 = np.zeros((rows,cols))
for i in range(rows):
    for j in range(cols):
        M3[i,j] = int(M2[i,j])
        
plt.matshow(M3)