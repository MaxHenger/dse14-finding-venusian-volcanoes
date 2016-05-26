# -*- coding: utf-8 -*-
"""
Created on Wed May 25 11:58:12 2016

@author: Yuyang
"""
import numpy as np
import sympy as sy
import matplotlib.pyplot as plt

"""definitions to get angle of attack"""
#put max acc conditions of a launcher into 1 block
def MaxAccConditions(launch_type, W_sc):
    if launch_type == 'block 1':
        #pointing direction
        """should change"""
        beta = 45.*np.pi/180 #np.linspace(90, 0, 91) #rad
        #angle of attack
        """should change"""
        alpha = 15.*np.pi/180 #rad
        #flight path angle
        """should change"""
        theta = 30.*np.pi/180 #rad
        
        #specific impulse
        #booster each
        """need verification"""
        Isp_booster = 286 #s
        #core stage
        Isp_core = 452.4 #s
        #thrust
        #core stage vacuum 100% throttle
        T_core =  232375.4*4/1.09 #kg
        #max acc thrust is when core stage drops
        T_maxacc = T_core #kg
        
        #gross lift off weight of launch vehicle
        #http://www.nasa.gov/pdf/664158main_sls_fs_master.pdf
        """this value needs verification"""
        W = 2650000. #2494758. #kg
        #glow of core stage
        #https://en.wikipedia.org/wiki/Space_Launch_System#cite_note-NASA_Fact-Sht-6
        """not secure site"""
        W_core = 979452 #kg
        #glow of booster each
        #http://www.spacelaunchreport.com/sls0.html
        """not secure site"""
        W_booster = 729800. #kg
        #max payload
        #SLS Mission Planners Guide (MPG) Overview, 02/2014
        W_pl = 20000. #kg
        #weight at max acceleration as soon as core stage drops, 6.2 km2/s2 c3
        W_maxacc = W - W_core - W_booster*2 - W_pl + W_sc
        
        #gravity at 2nd stage launch
        """not sure to use this value"""
        g = 9.6 #m/s2
        #density of air
        """should change with altitude"""
        rho = 1.225 #kg/m3, at 0 altitude
        
        #lift coefficient of launch vehicle, function of alpha
        Cl_lv = 2*np.pi*alpha
        #drag coefficent of launch vehicle
        """need to change"""
        Cd_lv = 0.5
        #velocity
        """need to change"""
        V = 50. #m/s
        #dynamic pressure 
        q = 0.5*rho*V**2 #Pa
        #aero area of launch vehicle
        """need to change"""
        S = 100 #m2
        #lift due to launch vehicle shape
        """not axial yet!"""
        L = Cl_lv*q*S #N
        #drag due to launch vehicle shape
        """not lateral yet!"""
        D = Cd_lv*q*S #N
        
        #max acceleration occurs as soon as core stage drops
        #SLS Mission Planners Guide (MPG) Overview, 02/2014
        a = 3.5*g
    return a, g, T_maxacc, W_maxacc, theta, D, L        

#get angle of attack from 
def LaunchAcceleration(launch_type):
    #import calcuation conditions
    a, g, T, W, theta, D, L = MaxAccConditions(launch_type, W_sc)[:]
    
    #axial acceleration
    """not usable for now"""
#    ax = g*(T/W - np.sin(theta + alpha) - D/W*np.cos(alpha) + L/W*np.sin(alpha))
    #simplified version, no aero
    theta = 0. #rad
    
    #solve this polynomial of alpha
    alpha = sy.symbols('alpha')
    ax_simp = a*sy.sin(alpha)
    alpha = sy.solve(g*(T/W - sy.sin(alpha)) - ax_simp, alpha)[0] #rad
    
    #get ax_simp value
    ax_simp = a*sy.sin(alpha)
    
    #lateral acceleration
    """not usable for now"""
#    az = g*(np.cos(theta + alpha) + D/W*np.sin(alpha) + L/W*np.cos(alpha))
    #simplified version, no aero
    az_simp = g*(sy.cos(theta + alpha)) #m/s2
    return ax_simp, alpha, az_simp, T, W

"""execution to get angle of attack"""
#choose launcher type
launch_type = 'block 1'
#real payload mass
W_sc = 2000. #kg
#weight of a wing
"""should get this from WingboxBoom"""
W_wing = 100. #kg
#print LaunchAcceleration(launch_type)[:]

"""Stress cases"""
#for fuselage
def FuselageLaunchLoad(L_fs):
    #import acc, alpha, T, W
    ax, alpha, az, T, W = LaunchAcceleration(launch_type)[:]
    
    #Weight at 2nd stage in aixal (x) and lateral (z) dir
#    Wx = W*np.sin(alpha) #kg
#    Wz = W*np.cos(alpha) #kg
    
    #resultant force at 3.5g
    Fx = W_pl*ax #N
    Fz = W_pl*az #N
    
    #moment varying with length of fuselage
    #start from middle to right, then mirror for left
    arm_aft = np.linspace(0, L_fs/2., 20)
    Mz_fs_aft = 0.5*Fz*arm_aft
    
    Mz = np.zeros(np.size(Mz_fs_aft) - 1)
    for i in range(np.size(Mz_fs_aft) - 1):
        Mz[i] = Mz_fs_aft[-1 - i]
    Mz = np.concatenate((Mz, Mz_fs_aft))
    return Fx, Fz, Mz, Mz_fs_aft
    
def FuselageLaunchStress(R, t_fs):
    Fx, Fz, Mz, Mz_fs_aft = FuselageLaunchLoad(L_fs)[:]
    
    Ixx = Iyy = np.pi*R**3*t_fs
    ymax = R
    
    sig_z_fs = Mz*ymax/Ixx
    sig_x_fs = Fx/np.pi*(R**2 - (R - t_fs)**2) #Pa, a number
    
    #plot the sigma
    plt.plot(np.linspace(0, L_fs, np.size(sig_z_fs)), sig_z_fs)
    plt.show()
    return sig_z_fs, sig_x_fs

"""need to optimise t_fs for different y location"""

"""exeution to get Mz"""
L_fs = 5.
W_pl = 1000.
R = 0.6
t_fs = 0.005
sig_z_fs, sig_x_fs = FuselageLaunchStress(R, t_fs)[:]
print FuselageLaunchStress(R, t_fs)[1]


#for wing, only mass of wing suffers 3.5g
def WingLaunchStress(W_wing):
    #import acc, alpha, T, W
    ax, alpha, az, T, W = LaunchAcceleration(launch_type)[:]
    
    #this axial force should be supported at root and tip
    Fx = W_wing*np.sin(alpha)*ax
    Fx_tip = Fx_root = Fx/2
    
    #same for lateral force dur to weight of wing at 3.5g
    Fz = W_wing*np.cos(alpha)*a
    Fz_tip = Fz_root = Fz/2
    
    return Fx_tip, Fx_root, Fz_tip, Fz_root
    
    
    
    
    
    
