# -*- coding: utf-8 -*-
"""
Created on Tue May 10 18:46:59 2016

This file contains the three entry methods( balistic, gliding and skipping) in
form of three classes. Each class has seperate methods for retrival of various
unique and general information. All contain a method for expected landing as
well as error bars indicating possible landing ranges. 

Ballistic class:

Gliding class:

Skipping class:

In general all functions accept and return values in SI units. An exception is
angles, these are entered in degrees. Also attempt to call the functions with
arrays as much as possible to reduce unneccesary calculations.

The general use of this file is as following:
    
    import Entry
    import Orbit

    orb = Orbit.Orbit('preliminary')
    orb.create(SemiMajor,Eccentricity,Inclination,AscentionAngle)
    orb.place(Location,Time)
    
    bal = Entry.Ballistic(orb)

@author: Julius
"""

# Given a certain orbit, entry method, landing site: calculate time to release
# orbit and release location, decent profile, decent acceleration, decent heat 
# flux, decent heat

#
#
#

import matplotlib.pyplot as plt
import Orbit
import Gravity
import numpy as np
import Utility as util
import Atmosphere
import control

atm = Atmosphere.Atmosphere("preliminary")



class Entry2D:
    """2D simulation to verify the 3D """
    def __init__(self,Orbit,mass,gravity=None):
        self.orbit=Orbit
        if type(gravity)==type(None):
            self.grav=self.orbit.grav
        else:
            self.grav=gravity

class Entry3D:
    """Ballistic entry profile: takes orbit as inital argument """
    def __init__(self,Orbit,mass,gravity=None):
        self.orbit=Orbit
        
    def __simulate__(self,dt):
        pass
        
        
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

## testing
import EntryCoef as ec
grav=Gravity.Gravity(radians=True)    
coeff = ec.coeff()

K=348.8 kg/m2
V0 = 11.2*1000 # m/s inital velocity 
y0 = 20*np.pi/180. # rad inital gamma, flight path angle, must be positive for convergence
X0 = 70.75 *np.pi/180. # rad inital heading
R0 = 300*1000+6052*1000 # m inital altitude
rho0 =util.scale_height(R0-6052*1000)[2]
tau0 = -106.58 *np.pi/180 # rad inital latitude
delta0= -22.3* np.pi/180. # rad inital longitude

ydot0 = -0.1*np.pi/180. # rad/s inital change in flight path angle
p0    = 0*np.pi/180.   # rad/s inital roll rate
q0    = 0*np.pi/180.    # rad/s inital pitch rate
r0    = 0*np.pi/180.    # rad/s inital yaw rate
a0    = 40*np.pi/180.    # rad inital pitch
b0    = 0*np.pi/180.    # rad inital yaw
m0    = 0*np.pi/180.    # rad inital roll

CL=0.1 #
L0=100  # N inital Lift
D0=100000  # N inital Drag
g0=grav(R0,tau0,delta0)

Ixx=10**-2
Iyy=10**-1
Izz=10**-2
Ixz=0.1

CDa= 0.1  # positive change in drag due to angle of attack
CDM= 0.1  # positive change in drag due to Mach
CLM= 0.1  # positive change in lift due to Mach
CLa= 0.2 # positive change in lift due to angle of attack
Clb= 0.1  # ?? change in roll due to sideslip
CmM=-0.1  # negative? change in pitch moment due to Mach 
Cma=-0.8  # negative change in pitch moment due to angle of attack
Cnb= 0.1  # positve change in yaw moment due to slide slip
CSb= 0  # positive? change in side force due to side slip
# Input stab
Clda= 0.1 # lift due to angle of attack
Cmda=-0.1 # pitch moment due to angle of attack
Cnda= 0 # Yaw moment due to angle of attack
Cndr= 0 # Yaw moment due to yaw rate

m=26029. # kg mass
S=110. # m2 surface area
b=13.
c=23

coeff.initial(V0,y0,ydot0,R0,p0,q0,r0,a0,b0,m0,rho0)
coeff.initial_misc(CL,L0,D0,g0)
coeff.initial_MoI(Ixx,Iyy,Izz,Ixz)
coeff.initial_stab(CDa,CDM,CLM,CLa,Clb,CmM,Cma,Cnb,CSb,Clda,Cmda,Cnda,Cndr)
coeff.initial_size(m,S,b,c)

MA = coeff.matrixA()
MB = coeff.matrixB()
MC = np.eye(len(MA))
MD = np.zeros(MB.shape)        
    
ssS=control.ss(MA,MB,MC,MD)
T=np.arange(0,10000,.1)
u=np.zeros( (len(T),6 ) )
u[0][:10]=0
u[1][1000:1100]=0

yout,T,xout = control.lsim(ssS,u,T,[V0,y0,R0,p0,q0,r0,a0,b0,m0])
y = yout.T
#plt.plot(T,y[0],label="V")
#plt.plot(T,y[1],label="y")
#plt.plot(T,y[2],label="R")
#plt.plot(T,y[3],label="p")
#plt.plot(T,y[4],label="q")
#plt.plot(T,y[5],label="r")
#plt.plot(T,y[6],label="a")
#plt.plot(T,y[7],label="b")
#plt.plot(T,y[8],label="m")
#plt.grid(True)
#plt.legend(loc=2)

figd=plt.figure("Entry Simulation Results", figsize=(20,11))
ncols=3
nrows=3
n=0
ax=[]
for i in range(nrows):
        for j in range(ncols):
            n+=1
            ax.append( figd.add_subplot(nrows*100+ncols*10+n) )
    
n=0
plt.rcParams.update({"font.size":14}) 
ax[n].plot(T, y[n], label = r"Velocity", color="blue", linestyle="-")
#ax[n].plot(T, roll_sim,label = r"Numerical improved", color="red", linestyle= "--" )
#ax[n].plot(T, roll_sim_old,label = r"Numerical original", color="green",linestyle=":", linewidth=2.5)
ax[n].set_title(r"Velocity")
ax[n].set_ylabel(r"Velocity in $m/s$", fontsize=16)

n=1
ax[n].plot(T, y[n], label = r"flight path", color="blue", linestyle="-")
#ax[n].plot(T, yaw_sim,label = r"Numerical improved", color="red", linestyle= "--" )
#ax[n].plot(T, yaw_sim_old,label = r"Numerical original", color="green",linestyle=":", linewidth=2.5)
ax[n].set_title(r"Flight path angle")
ax[n].set_ylabel(r"$\gamma$ in $rad$", fontsize=16)
n=2
ax[n].plot(T, y[n], label = r"Radius", color="blue", linestyle="-")
#ax[n].plot(T, rollrate_sim,label = r"Numerical improved", color="red", linestyle= "--" )
#ax[n].plot(T, rollrate_sim_old,label = r"Numerical original", color="green",linestyle=":", linewidth=2.5)
ax[n].set_title(r"Orbit Radius")
ax[n].set_ylabel(r"$R$ in $m$", fontsize=16)
n=3
ax[n].plot(T, y[n], label = r"roll rate", color="blue", linestyle="-")
#ax[n].plot(T, yawrate_sim,label = r"Numerical improved", color="red", linestyle= "--" )
#ax[n].plot(T, yawrate_sim_old,label = r"Numerical original", color="green",linestyle=":", linewidth=2.5)
ax[n].set_title(r"Roll rate")
ax[n].set_ylabel(r"$p$ in $rad / s$", fontsize=16)
n=4
ax[n].plot(T, y[n], label = r"pitch rate")
#ax[n].plot(T, aileron_input,label = r"Aileron", linestyle= "--" )
ax[n].set_title(r"pitch rate")
ax[n].set_ylabel(r"$q$ in $rad /s$", fontsize=16)
n=5
ax[n].plot(T, y[n], label = r"yaw rate")
#ax[n].plot(T, aileron_input,label = r"Aileron", linestyle= "--" )
ax[n].set_title(r"yaw rate")
ax[n].set_ylabel(r"$r$ in $rad / s$", fontsize=16)
n=6
ax[n].plot(T, y[n], label = r"angle of attack")
#ax[n].plot(T, aileron_input,label = r"Aileron", linestyle= "--" )
ax[n].set_title(r"angle of attack")
ax[n].set_ylabel(r"$\alpha$ in $rad$", fontsize=16)
n=7
ax[n].plot(T, y[n], label = r"yaw angle")
#ax[n].plot(T, aileron_input,label = r"Aileron", linestyle= "--" )
ax[n].set_title(r"yaw angle")
ax[n].set_ylabel(r"$\beta$ in $rad$", fontsize=16)
n=8
ax[n].plot(T, y[n], label = r"roll angle")
#ax[n].plot(T, aileron_input,label = r"Aileron", linestyle= "--" )
ax[n].set_title(r"roll angle")
ax[n].set_ylabel(r"$\mu$ in $rad$", fontsize=16)

    
for i,axd in enumerate(ax):
    axd.grid(True)
    axd.set_xlabel(r"Time in $s$", fontsize=14)
    if i==1:
        axd.legend(loc=3,prop={'size':12})
    elif i==4:
        axd.legend(loc=2,prop={'size':16})
    else:
        axd.legend(loc=2,prop={'size':12})
    

plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=1.4)
plt.show()



if __name__!="__main__":
    grav=Gravity.Gravity()    
    
    SemiMajor       = 10439.7   # km 
    Eccentricity    = 0.39177   # -
    Inclination     = 85.5      # deg
    AscentionAngle  = 272.69    # deg
    TrueAnomaly     = 180       # deg : at periapsis
    ArgumentPeri    = 168.5     # deg
    orb = Orbit.Orbit(grav,SemiMajor,Eccentricity,Inclination,AscentionAngle,ArgumentPeri,TrueAnomaly)
    
    mass=100
    
    ent = Entry2D(orb,mass,grav)
    ent.simulate(0.01)
