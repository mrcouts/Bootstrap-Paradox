# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from math import pi as Pi
from txt2py import *

A = txt2py("aquisicao_circulo_PIDSMCt_lambda35_phi20_k200_6.txt")

#aquisicao_circulo_PIDSMCx_lambda70_phi6_k70_0_3: e_quad = 0.912252649215 | tau_quad = 0.238929403956 | s1_quad = 4.83123738371 | s2_quad = 4.92973723999 | ex_quad = 0.693586588682 | ey_quad = 0.592572814091 | tau1_quad = 0.163618207999 | tau2_quad = 0.174115886943
#aquisicao_SMCt_lambda60_phi20_k363_3:            e_quad = 1.0865600936   | tau_quad = 0.237820339378 | s1_quad = 55.7256386861 | s2_quad = 52.1114451708 | ex_quad = 0.89867749944  | ey_quad = 0.610730373405 | tau1_quad = 0.162677093788 | tau2_quad = 0.173478174358
#aquisicao_circulo_SMCx_lambda60_phi6_k103_0_3:   e_quad = 1.27647821995  | tau_quad = 0.260967898131 | s1_quad = 8.35898181713 | s2_quad = 4.9935905567  | ex_quad = 0.971027682922 | ey_quad = 0.828554092984 | tau1_quad = 0.169215497632 | tau2_quad = 0.198671485664
#aquisicao_circulo_PIDx_lambda70_4:               e_quad = 1.28260640104  | tau_quad = 0.219666166256 | s1_quad = 6.09596142926 | s2_quad = 4.45005347346 | ex_quad = 1.14463260263  | ey_quad = 0.578701464488 | tau1_quad = 0.149533211028 | tau2_quad = 0.16091315483
#aquisicao_circulo_CTCt_lambda60_3:               e_quad = 1.61909535266  | tau_quad = 0.202338133235 | s1_quad = 45.0666924646 | s2_quad = 53.5093564972 | ex_quad = 1.31948403439  | ey_quad = 0.938313190784 | tau1_quad = 0.132682753043 | tau2_quad = 0.152761275217
#aquisicao_circulo_CTCx_lambda60_1:               e_quad = 1.63913254193  | tau_quad = 0.212001438004 | s1_quad = 5.99433632576 | s2_quad = 4.31327608202 | ex_quad = 1.358628447    | ey_quad = 0.91699740076  | tau1_quad = 0.140617598848 | tau2_quad = 0.158654658331
#aquisicao_circulo_PIDSMCt_lambda35_phi20_k200_6: e_quad = 3.34857383374  | tau_quad = 0.198843595669 | s1_quad = 38.0594077422 | s2_quad = 36.1432633763 | ex_quad = 1.6688967167   | ey_quad = 2.90305536788  | tau1_quad = 0.135978269066 | tau2_quad = 0.14508165246
#aquisicao_circulo_PIDt_lambda35_11:              e_quad = 4.45725019704  | tau_quad = 0.186064672078 | s1_quad = 34.1700477306 | s2_quad = 32.2174188766 | ex_quad = 2.3747188143   | ey_quad = 3.77197426714  | tau1_quad = 0.120771293755 | tau2_quad = 0.141542773748

#aquisicao_triangulo_SMCx_lambda60_phi6_k103_0_2:    e_quad = 1.930624, i_quad = 4.677893, ex_quad = 1.499759, s1_quad = 8.5627540, s2_quad = 15.5908030, ey_quad = 1.215743,  i1_quad = 3.047811,  i2_quad = 3.548737  
#aquisicao_triangulo_SMCt_lambda60_phi20_k363_2:     e_quad = 2.315927, i_quad = 4.790667, ex_quad = 2.073283, s1_quad = 90.622643, s2_quad = 96.8548660, ey_quad = 1.031995,  i1_quad = 3.142629,  i2_quad = 3.615850 
#aquisicao_triangulo_PIDSMCx_lambda60_phi6_k103_0_3: e_quad = 2.534649, i_quad = 4.721774, ex_quad = 0.985656, s1_quad = 9.2306470, s2_quad = 23.4585950, ey_quad = 2.335151,  i1_quad = 2.900244,  i2_quad = 3.726088  
#aquisicao_triangulo_CTCx_lambda60_4:                e_quad = 2.593759, i_quad = 4.709301, ex_quad = 2.198307, s1_quad = 10.335380, s2_quad = 15.6801630, ey_quad = 1.376601,  i1_quad = 3.014247,  i2_quad = 3.618262 
#aquisicao_triangulo_CTCt_lambda60_5:                e_quad = 2.607846, i_quad = 4.541315, ex_quad = 2.244732, s1_quad = 85.878242, s2_quad = 106.841187, ey_quad = 1.327417,  i1_quad = 2.874919,  i2_quad = 3.515448
#aquisicao_triangulo_PIDx_lambda70_6:                e_quad = 2.803854, i_quad = 4.803226, ex_quad = 1.434845, s1_quad = 10.004310, s2_quad = 23.8964730, ey_quad = 2.408903,  i1_quad = 2.968096,  i2_quad = 3.776425 
#aquisicao_triangulo_PIDSMCt_lambda25_phi90_k200_2:  e_quad = 7.304744, i_quad = 5.331062, ex_quad = 4.019952, s1_quad = 81.700752, s2_quad = 70.6879430, ey_quad = 6.099121,  i1_quad = 3.646859,  i2_quad = 3.888527 
#aquisicao_triangulo_PIDt_lambda25_2:                e_quad = 7.615739, i_quad = 4.938209, ex_quad = 3.926870, s1_quad = 70.306961, s2_quad = 66.9731520, ey_quad = 6.525272,  i1_quad = 3.150745,  i2_quad = 3.802461 


t_np = np.array([A[i][0] for i in range(len(A))])
x_np = np.array([A[i][1] for i in range(len(A))])
y_np = np.array([A[i][2] for i in range(len(A))])
xref_np= np.array([A[i][3] for i in range(len(A))])
yref_np= np.array([A[i][4] for i in range(len(A))])
ex_np   = np.array([A[i][5] for i in range(len(A))])
ey_np   = np.array([A[i][6] for i in range(len(A))])
i1_np   = np.array([A[i][7] for i in range(len(A))])
i2_np   = np.array([A[i][8] for i in range(len(A))])
u1_np   = np.array([A[i][9] for i in range(len(A))])
u2_np   = np.array([A[i][10] for i in range(len(A))])
s1_np   = np.array([A[i][11] for i in range(len(A))])
s2_np   = np.array([A[i][12] for i in range(len(A))])

t_np = 0.001*t_np
tau1_np = 0.055984*i1_np
tau2_np = 0.0566596*i2_np

"""T = 1000
Ta = 3

n2 = (2-1)*T/3
n8 = (8-1)*T/3

print n2
print n8

ex_quad = 0.0
ey_quad = 0.0
tau1_quad = 0.0
tau2_quad = 0.0
s1_quad = 0.0
s2_quad = 0.0

for i in range(n2,n8):
	ex_quad += ex_np[i]**2
	ey_quad += ey_np[i]**2
	tau1_quad += tau1_np[i]**2
	tau2_quad += tau2_np[i]**2
	s1_quad += s1_np[i]**2
	s2_quad += s2_np[i]**2

ex_quad = 1000*(2*ex_quad/(n8-n2))**0.5
ey_quad = 1000*(2*ey_quad/(n8-n2))**0.5
e_quad = (ex_quad**2 + ey_quad**2)**0.5
tau1_quad = (2*tau1_quad/(n8-n2))**0.5
tau2_quad = (2*tau2_quad/(n8-n2))**0.5
tau_quad = (tau1_quad**2 + tau2_quad**2)**0.5
s1_quad = (2*s1_quad/(n8-n2))**0.5
s2_quad = (2*s2_quad/(n8-n2))**0.5

print "e_quad =", e_quad, "| tau_quad =", tau_quad, "| s1_quad =", s1_quad, "| s2_quad =", s2_quad, "| ex_quad =", ex_quad, "| ey_quad =", ey_quad, "| tau1_quad =", tau1_quad, "| tau2_quad =", tau2_quad
"""

fig, ax = plt.subplots()
ax.plot(xref_np, yref_np, 'b', linewidth=1, label= 'Refer' + u'ê' 'ncia')
ax.plot(x_np,    y_np,    'r', linewidth=1, label='Trajet' + u'ó' + 'ria real')
plt.xlabel(r'$x[m]$')
plt.ylabel(r'$y[m]$')
plt.title('Trajet' + u'ó' + 'ria realizada')
ax.legend(loc=4, ncol=1, prop={'size': 10})
plt.savefig('xy.png')

plt.figure()
plt.plot(t_np, ex_np, 'r', linewidth=2)
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$e_x[m]$')
plt.title('Erro de posi' + u'ç' + u'ã' + 'o em fun'  + u'ç' + u'ã' 'o do tempo (coordenada x)')
plt.savefig('ex.png')

plt.figure()
plt.plot(t_np, ey_np, 'r', linewidth=2)
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$e_y[m]$')
plt.title('Erro de posi' + u'ç' + u'ã' + 'o em fun'  + u'ç' + u'ã' 'o do tempo (coordenada y)')
plt.savefig('ey.png')

plt.figure()
plt.plot(t_np, tau1_np, 'r', linewidth=2)
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$\tau_1[Nm]$')
plt.title('Torque aplicado ' + 'em fun'  + u'ç' + u'ã' 'o do tempo (atuador 1)')
plt.savefig('tau1.png')

plt.figure()
plt.plot(t_np, tau2_np, 'r', linewidth=2)
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$\tau_2[Nm]$')
plt.title('Torque aplicado ' + 'em fun'  + u'ç' + u'ã' 'o do tempo (atuador 2)')
plt.savefig('tau2.png')

plt.figure()
plt.plot(t_np, u1_np, 'r', linewidth=2)
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$u_1[V]$')
plt.title('Tens' + u'ã' + 'o aplicada ' + 'em fun'  + u'ç' + u'ã' 'o do tempo (atuador 1)')
plt.savefig('u1.png')

plt.figure()
plt.plot(t_np, u2_np, 'r', linewidth=2)
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$u_2[V]$')
plt.title('Tens' + u'ã' + 'o aplicada ' + 'em fun'  + u'ç' + u'ã' 'o do tempo (atuador 2)')
plt.savefig('u2.png')

plt.figure()
plt.plot(t_np, s1_np, 'r', linewidth=2)
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$s_1[rad/s^2]$')
plt.title('Vari' + u'á' + 'vel de escorregamento ' + 'em fun'  + u'ç' + u'ã' 'o do tempo (atuador 1)')
plt.savefig('s1.png')

plt.figure()
plt.plot(t_np, s2_np, 'r', linewidth=2)
plt.xlabel(r'$t[s]$')
plt.ylabel(r'$s_2[rad/s^2]$')
plt.title('Vari' + u'á' + 'vel de escorregamento ' + 'em fun'  + u'ç' + u'ã' 'o do tempo (atuador 2)')
plt.savefig('s2.png')