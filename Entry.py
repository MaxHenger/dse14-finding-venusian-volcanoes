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
        self.mass=mass
        if type(gravity)==type(None):
            self.grav=self.orbit.grav
        else:
            self.grav=gravity
    def simulate(self,dt):
        import sympy as sp
        V = sp.Symbol("V")
        Vdot = sp.Symbol("Vdot")
        L = sp.Symbol("L") # 1/2 rho V**2 S CL
        D = sp.Symbol("D") # 1/2 rho V**2 S CD
        gamma = sp.Symbol("gamma")
        gammadot = sp.Symbol("gammadot")
        rho = sp.Symbol("rho")
        R = sp.Symbol("R")
        Rdot = sp.Symbol("Rdot")
        
        g = self.grav(0,0,0)
        m = self.mass
        S = np.pi*3**2
        CL=0.1
        CD=0.05
        
        eq1 = -m*Vdot - D - m*g*np.sin(gamma)
        eq2 = -m*V*gammadot + L - m*g*np.cos(gamma)*(1-V**2/(g*R))
        eq3 = -Rdot +V*np.sin(gamma)
        eq4 = -L + 0.5*rho*V**2*S*CL
        eq5 = -D + 0.5*rho*V**2*S*CD
        
        ####
        #State Space
        ####
        
        # state variables
        # V,Vdot, L, D, gamma, gammadot, R, Rdot
        
        # Y variables
        # 
        x = [V,Vdot,L,Ldot,D,Ddot,gamma,gammadot,R,Rdot]
        A = 0
        B = 0
        C = np.eye(4)
        D = np.zeros(B.shape)        
            
        control.ss(A,B,C,D)

class Entry3D:
    """Ballistic entry profile: takes orbit as inital argument """
    def __init__(self,Orbit,mass,gravity=None):
        self.orbit=Orbit
        
    def __simulate__(self,dt):
        pass
        
        
        

## testing
import EntryCoef as ec
grav=Gravity.Gravity()    
coeff = ec.coef2D()

coeff.inital(V0,y0,ydot0,R0,q0,r0,a0,m0)
coeff.inital_misc(CL,L0,D0,M0,g0)
coeff.inital_MoI(Ixx,Iyy,Izz)
coeff.initial_size(m,S,b,c)
MA = coeff.matrixA()
MB = coeff.matrixB()

MC = np.eye(len(x))
MD = np.zeros(MB.shape)        
    
ssS=control.ss(MA,MB,MC,MD)
T=np.arange(0,10000,0.1)
u=np.zeros(len(T))
yout,T,xout = control.lsim(ssS,u,T)
y = yout.T
plt.plot(T,y[0])    

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
