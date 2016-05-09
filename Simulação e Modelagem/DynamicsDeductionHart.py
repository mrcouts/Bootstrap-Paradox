from Denavit import *

#Matriz de parametros de Denavit-Hartemberg
fDH_ = lambda q_,l_,lg_: Matrix([
[0, pi/2, 0          , q_[0]+pi/2, 0, 0            , lg_[0]       , 'R'],
[0, pi/2, l_[0]+l_[1], q_[1]+pi  , 0, -l_[1]+lg_[1], 0            , 'R'],
[0, pi/2, 0          , q_[2]+pi  , 0, 0            , lg_[2]       , 'R'],
[0, pi/2, l_[2]+l_[3], q_[3]+pi  , 0, -l_[3]+lg_[3], 0            , 'R'],
[0, pi/2, 0          , q_[4]+pi  , 0, 0            , lg_[4]       , 'R'],
[0, 0   , l_[4]+l_[5], q_[5]     , 0, 0            , -l_[5]+lg_[5], 'R']
])

R = Serial('6R', '', 6, fDH_, Matrix([0,1,0]))