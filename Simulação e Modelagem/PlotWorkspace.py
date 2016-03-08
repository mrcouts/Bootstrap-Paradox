import matplotlib.pyplot as plt
import numpy as np
from InKinNum import *


#
lx=0.25 # Largura do retangulo
ly=0.28 # Altura do retangulo
dl=0.0005 # discretizacao
nx=int(lx/dl) # quadradinhos em x
ny=int(ly/dl) # quadradinhos em y
M=zeros(ny,nx) # Malha da discretizacao
epsx=0.0004  # valor maximo do modulo do determinante x para o ponto ser considerado singular
epst=0.00007 # valor maximo do modulo do determinante y para o ponto ser considerado singular

# circunferencia para tracar o disco cuja trajetoria sera realizada
r = 0.085
x0 = 0.0
y0 = 0.16


# Preenchimento da matriz M
for i in range(ny):
    for j in range(nx):
        x=j*dl+0.5*dl
        y=i*dl+0.5*dl

        xn_=Matrix([x,y])
        thetan_=ftheta_(xn_)
        
        if im(thetan_[0])==0 and im(thetan_[1])==0:
            if abs(det(fJx_(xn_,thetan_)))<epsx or abs(det(fJt_(xn_,thetan_)))<epst:
                M[i,j]=2 #Se M=2 entao o ponto o ponto eh singular
            else:
                M[i,j]=1 #Se M=1 entao o ponto eh nao singular e pertence a area de trabalho
                
        if (x - x0)**2 + (y - y0)**2 < r**2:
            M[i,j]=3 #Se M=3 entao o ponto pertence ao disco
        #Se M=0 entao nao pertence a area de trabalho

#Inversao das linhas da matriz para que a origem do sistema de coordenadas fique embaixo
M2=M.extract(range(ny-1,-1,-1),range(nx))
pprint(M2)

#plotando a matriz
M3 = np.zeros((ny,nx))
for i in range(ny):
    for j in range(nx):
        M3[i,j] = int(M2[i,j])
        
plt.matshow(M3)



