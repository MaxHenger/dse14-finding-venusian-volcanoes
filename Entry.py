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
import EntryCoef as ec
import EntryVehicle as ev

atm = Atmosphere.Atmosphere("preliminary")

class Entry:
    def __init__(self,orbit,vehicle,coeff=None,Grav=None):
        self.orbit=orbit
        self.vehicle=vehicle
        self.youts=[]
        if coeff==None:
            self.coeff=ec.coeff()
        if Grav==None:
            self.grav=Gravity.Gravity()
        else:
            self.grav=Grav
            
    def simulate(self,length=100,len_block=10,dt=0.1):
        self.T = np.arange(0,length+dt,dt)
        count_blocks=float(length)/len_block
        self.youts=[]
        for i in range(0,count_blocks+1):
            T0=i*len_block
            u = np.zeros( (len(len_block),6 ) )
            init=self.youts[-1][-1].T # [V0,y0,R0,p0,q0,r0,a0,b0,m0]
            self.youts.append(self.simulate_block(u,init,T0,len_block,dt))
        
    def simulate_block(self,u,init,T0,len_block,dt):
        T = np.arange(T0,T0+len_block+dt,dt)
        MA = coeff.matrixA()
        MB = coeff.matrixB()
        MC = np.eye(len(MA))
        MD = np.zeros(MB.shape)        
        ssS=control.ss(MA,MB,MC,MD)
        yout,T,xout = control.lsim(ssS,u,T,init)        
        
    def _update_Vehicle(self):
        self.Vehicle.update()
        
    def _inital(self,V0,y0,ydot0,R0,p0,q0,r0,a0,b0,m0,rho0):
        self.coeff.initial(V0,y0,ydot0,R0,p0,q0,r0,a0,b0,m0,rho0)
    
    def _inital_misc(self,CL,L0,D0,g0):
        self.coeff.initial_misc(CL,L0,D0,g0)
    def _initial_MoI(self,Ixx,Iyy,Izz,Ixz):
        self.coeff.initial_MoI(Ixx,Iyy,Izz,Ixz)
    def _initial_stab(self,CDa,CDM,CLM,CLa,Clb,CmM,Cma,Cnb,CSb,Clda,Cmda,Cnda,Cndr):
        self.coeff.initial_stab(CDa,CDM,CLM,CLa,Clb,CmM,Cma,Cnb,CSb,Clda,Cmda,Cnda,Cndr)
    def _initial_size(self,m,S,b,c):
        self.coeff.initial_size(m,S,b,c)
        
    def plot_entry_raw(self,T=None,yout=None):
        if yout==None:
            yout=np.hstack(self.youts)
            T = self.T
        else: 
            y = yout.T
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

## testing



if __name__=="__main__":
    grav=Gravity.Gravity(radians=True)    
    
    SemiMajor       = 10439.7   # km 
    Eccentricity    = 0.39177   # -
    Inclination     = 85.5      # deg
    AscentionAngle  = 272.69    # deg
    TrueAnomaly     = 180       # deg : at periapsis
    ArgumentPeri    = 168.5     # deg
    orb = Orbit.Orbit(grav,SemiMajor,Eccentricity,Inclination,AscentionAngle,ArgumentPeri,TrueAnomaly)
    # get 
    # flightpath angel
    # flightpath rotation
    # Velocity
    # Radius
    # heading
    # lat
    # long
    # g
    # rho
    #
    coeff = ec.coeff()
    vehicle = ev.Test_vehicle()
    entrySim = Entry(orb,vehicle,coeff)
    K_ref=348.8 #kg/m2
    
    ### Position stuff
    V0 = 11.2*1000 # m/s inital velocity 
    R0 = 200*1000+6052*1000 # m inital altitude
    X0 = 70.75 *np.pi/180. # rad inital heading
    tau0 = -106.58 *np.pi/180 # rad inital latitude
    delta0= -22.3* np.pi/180. # rad inital longitude
    
    ### Misc stuff
    g0=grav(R0,tau0,delta0)
    rho0 =util.scale_height(R0-6052*1000)[2]
    
    ### angle stuff
    y0 = -10*np.pi/180. # rad inital gamma, flight path angle, must be positive for convergence
    ydot0 = -0.001*np.pi/180. # rad/s inital change in flight path angle
    
    m0    = 10*np.pi/180.    # rad inital roll
    p0    = 0.001*np.pi/180.   # rad/s inital roll rate
    
    a0    = 40*np.pi/180.    # rad inital pitch
    q0    = 0*np.pi/180.    # rad/s inital pitch rate
    
    b0    = 10*np.pi/180.    # rad inital yaw
    r0    = 0.001*np.pi/180.    # rad/s inital yaw rate
    
    
    #entrySim.plot_entry_raw()

    def func(x, y):
        return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2
    from scipy.interpolate import griddata
    grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]
    points = np.random.rand(1000, 2)
    values = func(points[:,0], points[:,1])
    grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')
    plt.imshow(grid_z2)
#    ### Aerodynamic stuff
#    L_D=0.37
#    CD=1.25#(K_ref/m*S)**-1
#    CL=0.45#L_D*CD #
#    
#    ### Inertia stuff
#    Ixx=1000
#    Iyy=1000
#    Izz=1000
#    Ixz=100
#    
#    
#    
#    ### Stability stuff
#    CDa= 1./45  # positive change in drag due to angle of attack
#    CLa= 1./45 # positive change in lift due to angle of attack
#    Cma= -0.0299  # negative change in pitch moment due to angle of attack
#    
#    CDM= -0.003  # ?positive change in drag due to Mach
#    CLM= -0.03  # positive change in lift due to Mach
#    CmM= 1e-100  # negative change in pitch moment due to Mach 
#    
#    Clb= -0.1  # Zero change in roll due to sideslip
#    Cnb= 0.0365  # positve change in yaw moment due to slide slip
#    CSb= -1  # ?positive change in side force due to side slip
#    # Input stab
#    Clda= 0 # lift due to angle of attack
#    Cmda= 0 # pitch moment due to angle of attack
#    Cnda= 0 # Yaw moment due to angle of attack
#    Cndr= 0 # Yaw moment due to yaw rate
