from Denavit import *

dof = 2

#Matriz de parametros de Denavit-Hartemberg
fDH_ = lambda q_,l_,lg_: Matrix([
[l_[0], 0, 0, q_[0], -l_[0]+lg_[0], 0, 0, 'R'],
[l_[1], 0, 0, q_[1], -l_[1]+lg_[1], 0, 0, 'R']
])

l_ = [0.1, 0.1]
lg_= [0.05, 0.05]
m_ = [0.1, 0.1]
I__= [Matrix([[0,0,0],[0,0.0001,0],[0,0,0.0001]]),Matrix([[0,0,0],[0,0.0001,0],[0,0,0.0001]])]

R = Serial('RR', '', dof, fDH_, Matrix([0,-1,0]), l_=l_, lg_=lg_, m_=m_, I__=I__)

subs0 = {R.q_.M_[0]:pi/6, R.q_.M_[1]:pi/4}
subs1 = {R.dq_.M_[0]:20.0, R.dq_.M_[1]:30.0}