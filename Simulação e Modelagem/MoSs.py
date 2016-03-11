# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 11:39:31 2016

@author: andre
"""
from MatrixSymbol import *

t = symbols('t')

class System(object):
    """
    Variaveis da classe:
    - _ : Dicionario
         -System Label
         -Subsystems Labels
         -Description
         -Naming Rules
         -Replacement Rules
         -Reference Motion (para linearizacoes)
         -Constraint Order (q:Order no MoSsPack for Mathematica)
         -Generalized Variable Definition Order (q:Def:Order no MoSsPack for Mathematica)
         -Debug Mode
         -Timer
         -q[]
         -q+[]
         -*q[]
         -_q[]
         -*c[]
         -__ (_c no MoSsPack for Mathematica)
         -f
         -*f
         -*q? (Booleana)
         -A
         -B
         -C
         -R
         -S
         -C? (Booleana)
         -S? (Booleana)
         -Explicit EOM? (Booleana)
         -_dq
         -Subsystems (dicionario de subsistemas)
    Metodos da classe:
    - Acoplamento ( acoplar filhos )
    """
    def __init__(self):
        pass
    
    def Acoplamento(self):
        pass
    
class CorpoRigido3D(object):
    def __init__(self, label, GravitationalField = Matrix([0,0, -symbols('g')]), ExternalActiveTorque = zeros(3,1), ExternalActiveForce = zeros(3,1)):
        self._ = {}
        self._['System Label'] = label
        self._['Subsystems Labels'] = []
        self._['Description'] = 'Corpo rigido ' + str(label)
        self._['Naming Rules'] = {}
        self._['Replacement Rules'] = {}
        self._['Reference Motion'] = {}
        self._['Constraint Order'] = 1
        self._['Generalized Variable Definition Order'] = 1
        self._['Debug Mode'] = False
        self._['Timer'] = False
        
        vel_ = SMatrix( Matrix([ Function('vx_' + str(label))(t), 
                                 Function('vy_' + str(label))(t), 
                                 Function('vz_' + str(label))(t)]) )
                                              
        w_ = SMatrix( Matrix([ Function('wx_' + str(label))(t), 
                               Function('wy_' + str(label))(t), 
                               Function('wz_' + str(label))(t) ]) )
                                              
        self._['q'] = {'1':  vel_ + w_ }
        self._['q+'] = {'1': 0} #SMatrix
        self._['*q'] = {'1': 0} #SMatrix
        self._['*q?'] = False
        self._['_q'] = {}
        self._['*c'] = {}
        self._['__'] = {} #Concatenar dicionarios de _q e Replacement Rules
        
        m_ = symbols('m_'+str(label))
        I_ = diag(symbols('Jx_'+str(label)), symbols('Jy_'+str(label)), symbols('Jz_'+str(label)) )
        
        
        self._['f'] = ( SMatrix( Matrix([ -m_*vel_.M_.diff(t) + ExternalActiveForce ]), vel_.rowl_) 
                      + SMatrix( Matrix([ -I_*w_.M_.diff(t) - w_.M_.cross(I_*w_.M_) + ExternalActiveTorque ]), w_.rowl_) )
        self._['A'] = 0 #SMatrix
        self._['R'] = SMatrix(1, self._['q']['1'].rowl_, self._['q']['1'].rowl_)
        self._['B'] = self._['A']
        self._['C?'] = False
        self._['C'] = SMatrix(1, self._['q']['1'].rowl_, self._['q']['1'].rowl_)
        self._['S?'] = False
        self._['S'] = self._['C']
        self._['*f'] = self._['C'].T()*self._['f']
        self._['Explicit EOM?'] = False
        self._['_dq'] = {}
        self._['Subsystems'] = {}
    
    def Acoplamento(self):
        pass
