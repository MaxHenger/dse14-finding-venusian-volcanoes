# -*- coding: utf-8 -*-
"""
Created on Mon May 30 10:38:58 2016

@author: Julius
"""

import numpy as np
import matplotlib.pyplot as plt
import control
import math

    
def u_impulse(T,amplitude=1,delay=0,weight=None):
    """T is time array, amplitude is number, delay is T index""" 
    if weight is None:
        weight = [1,0]
    impulse = np.vstack((np.zeros(len(T)),np.zeros(len(T)))).T
    impulse[:][delay] = [weight[0]*amplitude*(np.pi/180),weight[1]*amplitude*(np.pi/180)]
    return impulse

def u_step(T,amplitude=1,delay=0,weight=[1,0]):
    de=np.hstack((np.zeros(delay) , (weight[0]*amplitude*(np.pi/180)*np.ones(len(T)-delay))))
    dt=np.hstack((np.zeros(delay) ,                 weight[1]*np.ones(len(T)-delay)))
    step = np.vstack( (de,dt) ).T
    return step
    
def u_block(T,length=20,amplitude=1,delay=0,weight=None):
    if weight is None:
        weight = [1,0]
    block = np.vstack((np.zeros(len(T)),np.zeros(len(T)))).T
    block[:][delay:delay+length] = [weight[0]*amplitude*(np.pi/180),weight[1]*amplitude*(np.pi/180)]
    return block
    
def plot_symmetric(c,yout,T):
    u,alpha,gamma,q=yout.T
    figi=plt.figure("Symmetric simulation Results", figsize=(20,10))
    ax1 = figi.add_subplot(231)
    ax2 = figi.add_subplot(232)
    ax3 = figi.add_subplot(233)
    ax4 = figi.add_subplot(234)
    ax5 = figi.add_subplot(235)
    ax6 = figi.add_subplot(236)
    
    ax1.plot(T, u, label=r"$\hat u$")
    ax1.set_title(r"$\hat u$ versus time")
    ax1.grid(True)
    ax1.legend(loc=5)
    
    ax2.plot(T, u*c.V0+c.V0*math.cos(c.gamma0), label = r"$Velocity$")
    ax2.set_title(r"Velocity in X axis versus time")
    ax2.grid(True)
    ax2.legend(loc=5)
    
    ax3.plot(T, (alpha)*180/np.pi, label = r"$\alpha$")
    ax3.set_title(r"$\alpha$ versus time")
    ax3.grid(True)
    ax3.legend(loc=5)  
    
    ax4.plot(T, gamma*180/np.pi, label = r"$\theta$")
    ax4.set_title(r"$theta$ versus time")
    ax4.grid(True)
    ax4.legend(loc=5) 
    
    ax5.plot(T, q*c.V0/c.c, label = r"$q$")
    ax5.set_title(r"$q$ versus time")
    ax5.grid(True)
    ax5.legend(loc=5) 
    
    dth=[0]+[(gamma[i+1]-gamma[i])/(T[i+1]-T[i]) for i in range(0,len(T)-1)]
    ax6.plot(T,q*c.V0/c.c,label = r"$q$")
    ax6.plot(T, dth, label = r"d$\theta$")
    ax6.set_title(r"gradiant $\theta$ and q versus time")
    ax6.grid(True)
    ax6.legend(loc=5) 
    
    plt.plot()

def plot_asymmetric(c,yout,T):
    beta,varphi,p,r=yout.T
    figs=plt.figure("Asymmetric simulation Results", figsize=(20,10))
    ax1 = figs.add_subplot(331)
    ax2 = figs.add_subplot(332)
    ax3 = figs.add_subplot(333)
    ax4 = figs.add_subplot(334)
    ax5 = figs.add_subplot(335)
    ax6 = figs.add_subplot(336)
    
    ax1.plot(T, beta, label=r"$\beta$")
    ax1.set_title(r"$\beta$ versus time")
    ax1.grid(True)
    ax1.legend(loc=5)
    
    ax2.plot(T, varphi, label = r"$\varphi$")
    ax2.set_title(r"$\varphi$ versus time")
    ax2.grid(True)
    ax2.legend(loc=5)
    
    ax3.plot(T, p*2*c.V0/c.b, label = r"$p$")
    ax3.set_title(r"$p$ versus time")
    ax3.grid(True)
    ax3.legend(loc=5)  
    
    ax4.plot(T, r*2*c.V0/c.b, label = r"$r$")
    ax4.set_title(r"$r$ versus time")
    ax4.grid(True)
    ax4.legend(loc=5) 

    ax5.plot(T, r*2*c.V0/c.b, label = r"$r$")
    ax5.plot(T, p*2*c.V0/c.b, label = r"$p$")
    ax5.set_title(r"$r$ and $p$ versus time")
    ax5.grid(True)
    ax5.legend(loc=5) 
    
    
    plt.plot()
    

if __name__=="__main__":
    import stability
    import AircraftStabilityCoeff


    co = AircraftStabilityCoeff.coeff()
    canard,main,tail=stability.dummyWings()
    
    V0=100
    rho0=1.5940
    CL=main.cl
    CD=main.cd
    alpha0=1*np.pi/180.
    S=main.surface
    c=main.chord
    b=main.span
    e=main.oswald
    A=main.aspect
    xcg = 1.367
    xac = 1.23
    m=660
    Ixx=1000
    Iyy=5000
    Izz=1000
    propInc=1*np.pi/180
    propArm=4
    co._steady_conditions(V0,rho0,CD,CL,alpha0)
    co._aircraft_properties(b,c,A,S,e,m,Ixx,Iyy,Izz,xcg,xac,propInc,propArm)
    co._tail(tail)
    co._mainWing(main)
    co._canard(canard)
    co.deriv()
    
    CZu=0
    Cma=0
    Cmu=0
    co.delta_long(Cma=Cma,CZu=CZu,Cmu=Cmu)
    
    ssS=co.stateSpace(symmetric=True)
    T=np.arange(0,100,0.1)
    u=np.zeros((len(T),2))
    
    upgust = 0.#15. #m/s
    alpha = np.tan(upgust/V0)
    #u[1][0]=alpha
    
    frontgust = 10. #m/s
    du = frontgust/V0
    #u[0][0]=du
    # uhat, alpha, theta, q*c/V]
    init=[du,alpha,0,0]
    print init[1]*180/np.pi
    yout,T,xout=control.lsim(ssS,u,T,init)
    plot_impulse_symmetric_SS(co,yout,T)
