# -*- coding: utf-8 -*-
"""
Created on Tue May 31 14:08:51 2016

@author: Mathijs

Aircraft control surface sizing
"""

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Integration function for control surface effectiveness, based on ratio between control surface and lifting surface areas
def ControlSurfaceEffectiveness(ControlsurfaceToLiftingsurfaceRatio):
    ratio = [0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7]
    eff = [0,0.16,0.26,0.35,0.41,0.47,0.52,0.56,0.6,0.64,0.67,0.71,0.74,0.8]
    f = interp1d(ratio, eff)

    return f(ControlsurfaceToLiftingsurfaceRatio)


"""
Aileron Sizing

Inputs:
b = wingspan [m]
C_root = root chord wing [m]
S = wing area [m^2]
ins = inside aileron location in fraction of semi-span [-]
out = outside aileron location in fraction of semi-span [-]
C_L_alphaW = wing lift slope, with preferably the compressibility taken into account [1/rad]
eff = effectiveness of aileron [-] DEPENDS ON RATIO AILERON CHORD TO WING CHORD AND COULD BE ADDED AS FUNCTION OR FOUND IN TABLE 
V_stall = stall speed [m/s]
rho_stall = density at stall speed [kg/m^3]
S_vt = vertical tail area [m^2]
S_can = canard area [m^2]
S_ht = horizontal tail area [m^2]
I_xx = moment of inertia of aircraft [m^4]

Outputs:
time = time to roll to the desired bank angle of 45 degrees [s]
aileron characteristics [y_in,y_out,b_ail,C_ail,A_ail]
y_in = inside coordinate of aileron [m]
y_out = outside coordinate of aileron [m]
C_ail = aileron chord [m]
A_ail = aileron area [m]

"""

def AileronSizing(b,C_root,S,ins,out,C_L_alphaW,ailtochord,taper,maxdeflection,V_stall,rho_stall,S_vt,S_can,S_ht,I_xx):
    # Calculate distances for the inner edge and outer edge of the aileron
    y_in = ins*(b/2.)
    y_out = out*(b/2.)

    eff = ControlSurfaceEffectiveness(ailtochord)

    # Determine aircraft rolling moment coefficient slope with respect to aileron deflection
    C_l_da = ((2.*C_L_alphaW*eff*C_root)/(S*b))*(((y_out**2/2.)+((2./3.*(taper-1)/b))*y_out**3)-((y_in**2/2.)+((2./3.*(taper-1)/b))*y_in**3))

    # Determine aircraft rolling moment coefficient for maximum aileron deflection
    C_l = C_l_da*np.deg2rad(maxdeflection)

    # Assume minimum flight speed is 1.1 times the stall speed
    V_sizing = 1.3*V_stall

    # Determine aircraft rolling moment at the most critical flight situation (1.2*V_stall)
    L_A = 0.5*rho_stall*V_sizing**2*S*C_l*b
    
    # Resistance of wing sections to roll motion [Roughly average value from source]
    C_DR = 0.9
    
    # Drag center location [m]
    y_D = 0.4*(b/2.)

    # Determine steady state roll rate
    P_ss = np.sqrt((2.*L_A)/(rho_stall*(S+S_ht+S_vt+S_can)*C_DR*y_D**3))

    # Determine steady state bank angle
    BankAngle = (I_xx/(rho_stall*y_D**3*(S+S_ht+S_vt+S_can)*C_DR))*np.log(P_ss**2)
    
    # Determine roll acceleration
    RollRate = P_ss**2/(2.*BankAngle)

    # Determine time required to get to the desired roll angle of 45 degrees, and should be time<1.7 s
    time = np.sqrt((2.*np.deg2rad(30.))/RollRate)
    
    # Determine inner and outer aileron chord, taking assuming taper of aileron matches that of the wing
    C_ailIn = C_root*(1+(2.*((taper-1.)/b)*y_in))*ailtochord
    C_ailOut = C_root*(1+(2.*((taper-1.)/b)*y_out))*ailtochord
    
    # Determine aileron area
    A_ail = 2.*((C_ailIn+C_ailOut)/2.)*(y_out-y_in)
    
    return time,[y_in,y_out,C_ailIn,C_ailOut,A_ail]
    
#if __name__ == "__main__":
#    C_L_alphaW = 4.5 # 1/rad
#    b = 14.49
#    C_root = 1.604
#    ins = 0.61
#    out = 0.95
#    ailtochord = 0.40
#    taper = 0.8
#    S = 21.
#    maxdeflection = 20. # deg
#    
#    V_stall = 41.1556
#    rho_stall = 1.225
#    
#    S_vt = 4.2
#    S_can = 0.
#    S_ht = 5.3
#    
#    I_xx = 28000.
#
#    Time,AileronDim = AileronSizing(b,C_root,S,ins,out,C_L_alphaW,ailtochord,taper,maxdeflection,V_stall,rho_stall,S_vt,S_can,S_ht,I_xx)
#    print Time
#    print AileronDim    
    

'''
Elevator sizing program

Inputs:
elevtochord = elevator chord to horizontal tail chord [-]
C_L_ah = lift curve slope for horizontal tail [1/rad]
dynPressureRatio = ratio of velocities of horizontal tail to wing [-]
spanratio = ratio of elevator span over horizontal tail span [-]
arearatio = ratio of elevator area over horizotail tail area [-]
W = weight of aircraft [N]
S = Wing area [m^2]
chord = wing Mean Aerodynamic Chord [-]
C_Lalpha = lift curve slope [1/rad]
x_ac_w = wing aerodynamic center x-coordinate [m]
z_T = arm of thrust force with respect to c.g. (Positive if below c.g.) [m]

rho = density at altitude of interest [-]
V_array = array of flight speeds from stall speed to maximum speed at the altitude wanted[m/s]
T_max = max thrust at altitude [N]

C_L0 = aircraft lift coefficient for zero angle of attack [-]
C_m0 = aircraft moment coefficient for zero angle of attack [-]

x_cgAft = most aft c.g. location with respect to datum [m]
x_cgFW = most forward c.g. location with respect to datum [m]
dedalpha = downwash gradient for horizontal tail [-]
x_ac_h = horizontal tail aerodynamic center location [m]
A_h = horizontal tail area (both sides) [m^2]

Output:
plot of elevator deflection vs flight speed to check that maximum deflection is smaller than limit of +/- 20 deg
A_elevator = approximation of area of elevators [m^2]
The elevator ratios can also be used in the CAD for the control surface
'''

def ElevatorSizing(elevtochord,C_L_ah,dynPressureRatio,spanratio,arearatio,W,S,chord,C_Lalpha,x_ac_w,z_T,rho,\
                    V_array,T_max,C_L0,C_m0,x_cgAft,x_cgFW,dedalpha,x_ac_h,A_h):
    
    # Create array of thrust with same dimension as V_array
    T_array = np.ones(len(V_array))*T_max    
    
    # Determine elevator efficiency
    eff = ControlSurfaceEffectiveness(elevtochord)

    # Determine several coefficients 
    C_L_de = -C_L_ah*dynPressureRatio*arearatio*spanratio*eff   

    # Determine dynamic pressure for all velocities
    q_dyn = 0.5*rho*V_array**2 
    
    # Determine C_L1
    C_L1 = W/(q_dyn*S) 
    
    # Determine vertical tail volume for aft c.g.
    V_H_Aft = arearatio*((x_ac_h-x_cgAft)/chord)  

    # Determine C_m_de for the aft c.g.    
    C_m_deAft = -C_L_ah*dynPressureRatio*V_H_Aft*spanratio*eff

    # Determine vertical tail volume for forward c.g.
    V_H_FW = arearatio*((x_ac_h-x_cgFW)/chord) 
    
    # Determine C_m_de for the forward c.g.
    C_m_deFW = -C_L_ah*dynPressureRatio*V_H_FW*spanratio*eff   

    # Determine C_malpha for both the aft and forward c.g. position
    C_malphaAft = C_Lalpha*((x_cgAft-x_ac_w)/chord)-C_L_ah*(1-dedalpha)*(dynPressureRatio**2)*V_H_Aft
    C_malphaFW = C_Lalpha*((x_cgFW-x_ac_w)/chord)-C_L_ah*(1-dedalpha)*(dynPressureRatio**2)*V_H_FW  
    
    # Determine elevator deflections for the range of velocities
    de_Aft = ((((T_array*z_T)/(q_dyn*S*chord))+C_m0)*C_Lalpha+((C_L1-C_L0)*C_malphaAft))/((C_Lalpha*C_m_deAft)-(C_malphaAft*C_L_de))
    de_FW = ((((T_array*z_T)/(q_dyn*S*chord))+C_m0)*C_Lalpha+((C_L1-C_L0)*C_malphaFW))/((C_Lalpha*C_m_deFW)-(C_malphaFW*C_L_de))    

    # Create plot
    plt.plot(V_array,-np.rad2deg(de_Aft),'b.',markersize=8,label="Most aft cg")
    plt.plot(V_array,-np.rad2deg(de_FW),'gx',markersize=10,label="Most forward cg")
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.xlabel(r'Speed (m/s)',fontsize=16)
    plt.ylabel(r'$\delta_E (deg)$',fontsize=18)
    plt.legend()
    plt.grid()
    plt.show()    
    
    A_elevator = arearatio*A_h
    
    return A_elevator

# Elevator check
if __name__ == "__main__":
    
    elevtochord = 0.49
    dynPressureRatio = 0.96 
    C_L_ah = 4.3
    spanratio = 1.
    arearatio = 16./70.
    W = 20000.*9.81
    S = 70.
    chord = 2.96
    C_Lalpha = 5.2
    x_ac_w = 0.5
    z_T = -0.3
    rho = 1.225
    V_array = np.linspace(44.,185.,100)
    T_max = 56000.
    C_L0 = 0.24
    C_m0 = 0.05
    x_cgAft = 0.8
    x_cgFW = 0.2
    
    dedalpha = 0.454
    x_ac_h = 12.
    
    
    ElevatorSizing(elevtochord,C_L_ah,dynPressureRatio,spanratio,arearatio,W,S,chord,C_Lalpha,x_ac_w,z_T,rho,\
                    V_array,T_max,C_L0,C_m0,x_cgAft,x_cgFW,dedalpha,x_ac_h,A_h)

'''
Rudder sizing program

Inputs:
rudtochord = rudder chord to vertical tail chord ratio [-]
arearatio = rudder area to vertical tail area ratio [-]
spanratio = rudder span to vertical tail span ratio [-]
C_L_aV = lift curve slope for vertical tail airfoil [1/rad]
x_ax_V = vertical tail aerodynamic center location with respect to datum at nose [m]
x_cg = center of gravity location with respect to datum at nose [m] RECOMMENDED TO USE MOST AFT C.G. BECAUSE THEN MOMENT ARM IS SMALLEST
dynPressureRatio = ratio of dynamic pressures over the vertical tail with respect to main wing [-]
b = wing span [m]
V_stall = aircraft stall speed [m/s]
S = wing area [m^2]

Output:
delta_R = rudder deflection required to generate the moment to counteract engine failure [deg]
This value should be checked with the maximum deflection allowed, and if required, the rudder characteristics should be adjusted.
Once satisfied, the different ratios used for the input can be utilized in the CAD of the control surface
'''

def RudderSizing(rudtochord,arearatio,spanratio,C_L_aV,x_ac_V,x_cg,dynPressureRatio,b,V_stall,S,ThrustEngines,y_T,EnginesOperative):
    
    # Determine rudder efficiency
    eff = ControlSurfaceEffectiveness(rudtochord)    
    
    # Determine l_V
    l_V = x_ac_V-x_cg
    # Vertical tail volume coefficient
    V_V = (l_V/b)*arearatio    

    # Determine coefficient yaw moment with respect to rudder deflection
    C_n_dR = -C_L_aV*V_V*dynPressureRatio*eff*spanratio
    
    # Set the design speed for the rudder sizing
    V_design = 1.0*V_stall    
    
    # Determine dynamic pressure
    q_dyn = 0.5*rho*V_design**2
    
    delta_R = np.rad2deg(np.sum(ThrustEngines*EnginesOperative*y_T)/(-q_dyn*S*b*C_n_dR))    
    
    print "Rudder deflection (deg)"
    return delta_R
    
#if __name__ == "__main__":
#    rudtochord = 0.3
#    arearatio = 26./125.
#    spanratio = 1.
#    C_L_aV = 4.5
#    x_ac_V = 18.
#    x_cg = 0.
#    dynPressureRatio = 0.97
#    b = 34.
#    V_stall = 56.59
#    rho = 1.225
#    
#    # Thrust force created by engines (left to right when looking from aircraft tail to nose)    
#    ThrustEngines = np.array([116000.,116000.])
#    # Moment arms of engine thrust forces (left of c.g. positive when looking from aircraft tail to nose)
#    y_T = np.array([6.,-6.])
#    # Set engines operative or not, 1 is operative, 0 is inoperative
#    EnginesOperative = np.array([1.,0.])
#    
#    S = 125.
#    
#    print RudderSizing(rudtochord,arearatio,spanratio,C_L_aV,x_ac_V,x_cg,dynPressureRatio,b,V_stall,S,ThrustEngines,y_T,EnginesOperative)
