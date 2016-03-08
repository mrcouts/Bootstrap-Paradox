from InvKinDeduction import *


# Dicionario dos paramentros geometricos do mecanismo
par = {l0: 0.05, l1:0.12, l2:0.15}

###                                    ###
# Solucao Numerica da Cinematica Inversa #
###                                    ###

fE1=lambda xn_: E1.subs(par).subs({x_[i]: xn_[i] for i in range(2)}).evalf()
fF1=lambda xn_: F1.subs(par).subs({x_[i]: xn_[i] for i in range(2)}).evalf()
fG1=lambda xn_: G1.subs(par).subs({x_[i]: xn_[i] for i in range(2)}).evalf()

fE2=lambda xn_: E2.subs(par).subs({x_[i]: xn_[i] for i in range(2)}).evalf()
fF2=lambda xn_: F2.subs(par).subs({x_[i]: xn_[i] for i in range(2)}).evalf()
fG2=lambda xn_: G2.subs(par).subs({x_[i]: xn_[i] for i in range(2)}).evalf()

ftheta_1 = lambda xn_: 2*atan((-fF1(xn_)-sqrt(fE1(xn_)**2+fF1(xn_)**2-fG1(xn_)**2))/(fG1(xn_)-fE1(xn_))) if fG1(xn_)!=fE1(xn_) else 2*atan(-fE1(xn_)/fF1(xn_)) 
ftheta_2 = lambda xn_: 2*atan((-fF2(xn_)-sqrt(fE2(xn_)**2+fF2(xn_)**2-fG2(xn_)**2))/(fG2(xn_)-fE2(xn_))) if fG2(xn_)!=fE2(xn_) else 2*atan(-fE2(xn_)/fF2(xn_))
ftheta_  = lambda xn_: Matrix([ftheta_1(xn_),ftheta_2(xn_)])


#Caculo das posicoes
def fx_(thetan_):
    """###                                    ###
    # Solucao Numerica da Cinematica Direta  #
    ###                                    ###"""
    xc1= (a1_[0]).subs(par).subs({theta_[i]: thetan_[i] for i in range(2)}).evalf()
    yc1= (a1_[1]).subs(par).subs({theta_[i]: thetan_[i] for i in range(2)}).evalf()
    xc2= (a2_[0]).subs(par).subs({theta_[i]: thetan_[i] for i in range(2)}).evalf()
    yc2= (a2_[1]).subs(par).subs({theta_[i]: thetan_[i] for i in range(2)}).evalf()
    r = par[l2]

    xm=(xc1+xc2)/2.0
    ym=(yc1+yc2)/2.0
    dx=xc1-xc2
    dy=yc1-yc2
    
    lamb = sqrt((r**2)/(dx**2 + dy**2) - 0.25)
    
    xp = xm - lamb*dy
    yp = ym + lamb*dx
    
    return Matrix([xp,yp])
    

#Calculo dos angulos intermediarios
def fthetai_(xn_,thetan_):

    x1 = par[l0] + par[l1]*cos(thetan_[0])
    y1 = par[l1]*sin(thetan_[0])
    x2 = -par[l0] - par[l1]*cos(thetan_[1])
    y2 = par[l1]*sin(thetan_[1])
    thetai_1 =  (atan2(xn_[1] - y1, xn_[0] - x1) - thetan_[0]).evalf()
    thetai_2 = (Pi-atan2(xn_[1] - y2, xn_[0] - x2) - thetan_[1]).evalf()
    return Matrix([thetai_1, thetai_2])

#Cinematica inversa de todas as coordenadas generalizadas
def fInvKin(xn_):
    thetan_ = ftheta_(xn_)
    thetain_= fthetai_(xn_,thetan_)
    qn_=Matrix([xn_[0], xn_[1],thetan_[0], thetain_[0], thetan_[1], thetain_[1]])
    return qn_
    
#Cinematica direta de todas as coordenadas generalizadas
def fDirKin(thetan_):
    xn_ = fx_(thetan_)
    thetain_=fthetai_(xn_,thetan_)
    qn_=Matrix([xn_[0], xn_[1],thetan_[0], thetain_[0], thetan_[1], thetain_[1]])
    return qn_

###                             ###
# Calculo Numerico dos Jacobianos #
###                             ###

fJx_=lambda xn_, thetan_: Jx_.subs(par).subs({x_[i]: xn_[i] for i in range(2)}).subs({theta_[i]: thetan_[i] for i in range(2)}).evalf()
fJt_=lambda xn_, thetan_: Jt_.subs(par).subs({x_[i]: xn_[i] for i in range(2)}).subs({theta_[i]: thetan_[i] for i in range(2)}).evalf()
