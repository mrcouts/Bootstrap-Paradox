from sympy import *
from math import pi as Pi
init_printing(use_unicode = True)
t=symbols('t')

#parametros do mecanismo
l0,l1,l2=symbols('l0,l1,l2')

#Coordenadas generalizadas
theta_ = Matrix([Function('theta_1')(t),Function('theta_2')(t) ]) #Angulos dos motores
x_     = Matrix([Function('x')(t),Function('y')(t)]) #Coordenadas do efetuador
q_     = Matrix([x_,theta_]) #Vetor com todas as coordenadas generalizadas

####                                 ####
#    Deducacao da Cinematica Inversa    #
####                                 ####

##                             ##
#     Cinematica de Posicao     #
##                             ##

#Coordenadas dos elos intermediarios em funcao de theta
a1_= Matrix([ l0+l1*cos(theta_[0]),l1*sin(theta_[0])])
a2_= Matrix([-l0-l1*cos(theta_[1]),l1*sin(theta_[1])])

#Vinculos de posicao
Eq1=(x_[0]-a1_[0])**2+(x_[1]-a1_[1])**2-l2**2
Eq2=(x_[0]-a2_[0])**2+(x_[1]-a2_[1])**2-l2**2
phi_ = Matrix([Eq1,Eq2]) #Vetor com as equacoes vinculares

#Solucao analitica das equacoes vinculares

E1=phi_[0].expand().simplify().diff(cos(theta_[0]))
F1=phi_[0].expand().simplify().diff(sin(theta_[0]))
G1=(phi_[0]-E1*cos(theta_[0])-F1*sin(theta_[0])).simplify()

E2=phi_[1].expand().simplify().diff(cos(theta_[1]))
F2=phi_[1].expand().simplify().diff(sin(theta_[1]))
G2=(phi_[1]-E2*cos(theta_[1])-F2*sin(theta_[1])).simplify()

##                             ##
#    Cinematica de Velocidade   #
##                             ##

Jx_=phi_.jacobian(x_)
Jt_=phi_.jacobian(theta_)
