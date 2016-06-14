
import numpy as np
from matplotlib import pyplot as plt


#g0=8.87
ge=9.81
rho0=65.
H=15.9*1000
beta=1/H
GM = 0.32486*10**6*10**9

h       = np.arange(0,300000)#70000. # m 
m       = (2000+700)  # kg
CD      = 1.18  # - # 0.7 for ballistic
S       = 23.559#np.pi*(4.6/2)**2     # m26051
gamma   =-10 *np.pi/180 # rad
Ve      = 13000.
Re      = 6052.*1000

l_cone = 2
r_nose = 0.5
r_base = 4.6
angle = np.tan((r_base-r_nose)/l_cone) 
print angle*180/np.pi

CD = (1-np.sin(angle)**4)*(r_nose/r_base)**2+2*np.sin(angle)**2*(1-(r_nose/r_base)**2*np.cos(angle)**2)
print CD

def g(h):
    g = GM/(Re+h)**2
    return g

def ballistic():
    V = Ve*np.e**(0.5*g(h)*rho0*np.e**(-beta*h)/(  abs(m*g(h))/(CD*S) * beta * np.sin(gamma) ) )
    a = beta * np.sin(gamma) * Ve**2 * (V/Ve)**2 * np.log(V/Ve)
    a_der = (np.insert(V,len(V),V[-1] )-np.insert(V,0,V[0]) )[:len(V)]*V*gamma
    return a,a_der,V
    
def plot_ballistic():
    aout,a_der,V = ballistic()
    #afilt = aout[~np.isnan(aout)]
    #print max(afilt/ge)

    figb=plt.figure("Ballistic Entry", figsize=(18,10))
    ax2 = figb.add_subplot(221)
    ax3 = figb.add_subplot(223)
    ax4 = figb.add_subplot(222)
    ax5 = figb.add_subplot(224)
    
#    vely = np.sin(np.deg2rad(gamma))*V
#    T  = h/vely
#    T-=T[-1]
#    print len(vely)
#    print len(T)
#    print V
    
    ax2.plot(V,h, label = r"$Velocity$")
    ax2.set_title(r"Velocity versus height", fontsize=14)
    ax2.set_xlabel(r"Velocity in $\frac{m}{s}$", fontsize=14)
    ax2.set_ylabel(r"height in $m$", fontsize=14)
    ax2.grid(True)
    
#    ax4.plot(T,h, label = r"$alt$")
#    ax4.set_title(r"Time vs Height", fontsize=14)
#    ax4.set_xlabel(r"Time in $s$", fontsize=14)
#    ax4.set_ylabel(r"height in $m$", fontsize=14)
#    ax4.grid(True)
    
    ax3.plot(aout/ge,h, label = r"$acceleration$")
    ax3.plot(-a_der/ge,h, label = r"$acceleration der$") 
    ax3.set_title(r"Acceleration in g versus height", fontsize=14)
    ax3.set_xlabel(r"Acceleration a/g", fontsize=14)
    ax3.set_ylabel(r"height in $m$", fontsize=14)
    ax3.legend(loc=2)
    ax3.grid(True)
    
    figb.tight_layout()
    figb.show()


plot_ballistic()

