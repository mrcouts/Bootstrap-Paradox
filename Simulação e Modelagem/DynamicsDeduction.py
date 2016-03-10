from Denavit import *


dof = 2

class PontualBody2D(object):
    """ Modelo de um ponto no espaco """
    def __init__(self):
        self.q_ = SMatrix( Matrix([Function('x')(t),Function('y')(t)]) )
        self.dq_= SMatrix( self.q_.M_.diff(t) )
        self.C_ = SMatrix(1, self.dq_.rowl_, self.dq_.rowl_)
        self.M_ = SMatrix( zeros(2), self.dq_.M_, self.dq_.M_ )
        self.v_ = SMatrix( zeros(2,1), self.dq_.M_)
        self.g_ = SMatrix( zeros(2,1), self.dq_.M_)

#Matriz de parametros de Denavit-Hartemberg
fDH_ = lambda q_,l_,lg_: Matrix([
[l_[0], 0, 0, q_[0], -l_[0]+lg_[0], 0, 0, 'R'],
[l_[1], 0, 0, q_[1], -l_[1]+lg_[1], 0, 0, 'R']
])

RR1 = Serial('RR1', '1', dof, fDH_, Matrix([0,-1,0]))
RR2 = Serial('RR2', '2', dof, fDH_, Matrix([0,-1,0]))
P   = PontualBody2D()

## Acoplamento dos subsistemas ##

#coordenadas generalizadas
qh_ = P.q_
qo_ = RR1.q_ + RR2.q_
q_ = qh_ + qo_ 

#velocidades generalizadas
ph_ = P.dq_
po_ = RR1.p_ + RR2.p_
p_ = ph_ + po_

#
rhoh_ = P.dq_
rhoo_ = RR1.dq_ + RR2.dq_
rho_ = rhoh_ + rhoo_

H1 = H('y',0,symbols('l_0'),0,0)
H2 = H('y',pi,-symbols('l_0'),0,0)

#Vinculos de posicao:
_q_ = SMatrix( Matrix([ qh_.M_ - (H1*Matrix([RR1.o__[2], Matrix([1]) ]))[0:2,0], qh_.M_ - (H2*Matrix([RR2.o__[2], Matrix([1]) ]))[0:2,0] ]) , list((RR1.dq_ + RR2.dq_).M_ ))

Jh_ = SMatrix( _q_.M_.jacobian(qh_.M_), _q_.rowl_, list(qh_.M_.diff(symbols('t'))) )
Jo_ = SMatrix( _q_.M_.jacobian(qo_.M_), _q_.rowl_, list(qo_.M_.diff(symbols('t'))) )
J_ = Jh_ + Jo_

Cch_ = SMatrix(1, ph_.rowl_, ph_.rowl_ ) + (-1*Jo_.inv()*Jh_).simplify()
C_ = ( (P.C_ + RR1.C_ + RR2.C_)*Cch_ ).simplify()
dC_ = SMatrix( C_.M_.diff(symbols('t')), C_.rowl_, C_.coll_ )
A_ = J_ + RR1.A_ + RR2.A_
Ah_ = A_.extract(A_.rowl_, ph_.rowl_)
Ao_ = A_.extract(A_.rowl_, po_.rowl_)
b_ = (-1*J_.diff(t)*(P.dq_ + RR1.dq_ + RR2.dq_)).simplify() + RR1.b_ + RR2.b_

M_ = RR1.M_ + RR2.M_ + P.M_
v_ = RR1.v_ + RR2.v_ + P.v_
g_ = RR1.g_ + RR2.g_ + P.g_

uh_ = SMatrix( Matrix([symbols('tau_1'),symbols('tau_2')]), [RR1.dq_.rowl_[0],RR2.dq_.rowl_[0]] )
u_ = SMatrix( 0, rho_.rowl_ ) + uh_
u2_ = SMatrix( 0, p_.rowl_ ) + uh_
U_ = SMatrix(u2_.M_.jacobian(uh_.M_), u2_.rowl_, uh_.rowl_)

Cch_T_u_ = Cch_.T()*u_
Z_ = SMatrix( Cch_T_u_.M_.jacobian(uh_.M_), Cch_T_u_.rowl_, uh_.rowl_ )