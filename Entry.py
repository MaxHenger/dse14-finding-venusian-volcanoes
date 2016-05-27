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
grav=Gravity.Gravity()    
coeff = ec.coef2D()

V0 = 7435.5 # m/s inital velocity 
y0 = -1.43*np.pi/180. # rad inital gamma, flight path angle
X0 = 70.75 *np.pi/180. # rad inital heading
R0 = 300*1000 # m inital height
tau0 = -106.58 *np.pi/180 # rad inital latitude
delta0= -22.3* np.pi/180. # rad inital longitude

ydot0=-.001 # rad/s inital change in flight path angle
p0=0 # rad/s inital roll rate
q0=0 # rad/s inital pitch rate
r0=0 # rad/s inital yaw rate
a0=0 # rad inital pitch
b0=0 # rad inital yaw
m0=0 # rad inital roll

CL=0 #
L0=0  # N inital Lift
D0=0  # N inital Drag
g0=grav(R0,tau0,delta0)

Ixx=10**-4
Iyy=10**-3
Izz=10**-4
Ixz=0

CDa= 0.1  # positive change in drag due to angle of attack
CDM= 0.1  # positive change in drag due to Mach
CLM= 0.1  # positive change in lift due to Mach
CLa= 0.1  # positive change in lift due to angle of attack
Clb= 0.1  # ?? change in roll due to sideslip
CmM=-0.1  # negative? change in pitch moment due to Mach 
Cma=-0.8  # negative change in pitch moment due to angle of attack
Cnb= 0.1  # positve change in yaw moment due to slide slip
CSb= 0.1  # positive? change in side force due to side slip
# Input stab
Clda=0.1
Cmda=0.1
Cnda=0.1
Cndr=0.1

m=26029. # kg mass
S=110. # m2 surface area
b=13.
c=S/b

coeff.initial(V0,y0,ydot0,R0,p0,q0,r0,a0,b0,m0)
coeff.initial_misc(CL,L0,D0,g0)
coeff.initial_MoI(Ixx,Iyy,Izz,Ixz)
coeff.initial_stab(CDa,CDM,CLM,CLa,Clb,CmM,Cma,Cnb,CSb,Clda,Cmda,Cnda,Cndr)
coeff.initial_size(m,S,b,c)

MA = coeff.matrixA()
MB = coeff.matrixB()
MC = np.eye(len(MA))
MD = np.zeros(MB.shape)        
    
ssS=control.ss(MA,MB,MC,MD)
T=np.arange(0,10000,0.1)
u=np.zeros( (len(T),6 ) )
u[0][:1]=1
u[1][1000:1100]=1

yout,T,xout = control.lsim(ssS,u,T)
y = yout.T
plt.plot(T,y[0])
plt.plot(T,y[1])
plt.plot(T,y[2])
plt.plot(T,y[3])
plt.plot(T,y[4])
plt.plot(T,y[5])
plt.plot(T,y[6])
plt.plot(T,y[7])
plt.plot(T,y[8])    



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
