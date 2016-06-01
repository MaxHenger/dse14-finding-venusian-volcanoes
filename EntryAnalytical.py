
import numpy as np
from matplotlib import pyplot as plt


g0=8.87
rho0=65.
H=15.9*1000
beta=1/H
GM = 0.32486*10**6*10**9

h       = np.arange(0,300000)#70000. # m 
m       = 2000+700   # kg
CL      = 0.2
CD      = 0.7   # - # 0.7 for ballistic
S       = np.pi*(4.6/2)**2     # m26051
gamma   =-30 *np.pi/180 # rad
Ve      = 13000.
Re      = 6052.*1000

def g(h):
    g = GM/(Re+h)**2
    return g

def ballistic():
    V = Ve*np.e**(0.5*g(h)*rho0*np.e**(-beta*h)/(  abs(m*g(h))/(CD*S) * beta * np.sin(gamma) ) )
    a = beta * np.sin(gamma) * Ve**2 * (V/Ve)**2 * np.log(V/Ve)
    return a/9.81,V
    

def gliding():
    V = Ve * np.sqrt( (m*g(h)/S)/CL / ( 0.5*rho0*np.e**(-beta*h)*Ve**2 + (m*g(h)/S)/CL )  ) 
    Vc = np.sqrt(g(h)*(Re+h))
    #gamma = (-1/(beta)*2/(CL/CD)*g(h)/V**2)
    rho = 2* (m*g(h)/S)/CL*(1/V**2 - 1/Vc**2)
    gamma = np.arcsin( - 4/(beta*rho)*(m*g(h)/S)/CL * 1/V**3 * CD/CL*g(h)*(1-V**2/Vc**2)/V  )    
    #a = +CD/CL*(1-V**2/Vc**2 - 2/(beta*Re)*1/V**2/Vc**2) 
    a_der = (np.insert(V,len(V),V[-1] )-np.insert(V,0,V[0]) )[:len(V)]*V*gamma
    #a = -CD/(CL* 1/(1-V**2/Vc**2)) +(gamma)
    a = -CD/(CL* 1/(1-V**2/Vc**2))
    return a,a_der/9.81,gamma,Vc,V,rho
    
#
gout,V = ballistic()
gfilt = gout[~np.isnan(gout)]
print max(gfilt)

glide=False

if glide:
    a,a_der,gamma,Vc,V,rho = gliding()
    figs=plt.figure("Step_response", figsize=(18,10))
    ax2 = figs.add_subplot(221)
    ax3 = figs.add_subplot(223)
    ax4 = figs.add_subplot(222)
    ax5 = figs.add_subplot(224)
    
    ax2.plot(h,V, label = r"$Velocity$")
    ax2.set_title(r"Velocity versus height", fontsize=14)
    ax2.set_ylabel(r"Velocity in $\frac{m}{s}$", fontsize=14)
    ax2.set_xlabel(r"height in $s$", fontsize=14)
    ax2.grid(True)
    
    ax3.plot(h,a, label = r"$acceleration$")
    ax3.plot(h,gamma, label = r"$gamma$") 
    ax3.plot(h,a_der, label = r"$acceleration der$") 
    ax3.set_title(r"Acceleration in g versus height", fontsize=14)
    ax3.set_ylabel(r"Acceleration a/g", fontsize=14)
    ax3.set_xlabel(r"height in $s$", fontsize=14)
    ax3.legend(loc=2)
    ax3.grid(True)

