# -*- coding: utf-8 -*-
"""
Created on Wed May 25 17:04:34 2016

@author: Mathijs

Heat transfer functions

"""

import numpy as np

"""
Conduction function

Inputs:
k = thermal conductivity [W/(m*K)]
Area = area over which heat is transferred [m^2]
T_h = highest temperature [K or Celsius]
T_l = lowest temperature [SAME UNIT AS T_h!!!]
t = time [s]
d = thickness [m]

Output:
Q_cond = conducted heat over time t [J]
"""


def HeatTransferConduction(k,Area,T_h,T_l,t,d):
    
    Q_cond = (k*Area*(T_h-T_l)*t)/d
    return Q_cond    
    
"""
Convection function

Inputs:
H_c = heat transfer coefficient [W/(m^2*K)]
Area = area over which heat is transferred [m^2]
T_h = highest temperature [K or Celsius]
T_l = lowest temperature [SAME UNIT AS T_h!!!]
t = time [s]

Output:
Q_conv = convected heat over time t [J]
"""    
    
def HeatTransferConvection(H_c,Area,T_h,T_l,t):
    
    Q_conv = H_c*Area*(T_h-T_l)*t
    return Q_conv    
    
"""
Radiation functions

Inputs:
T_h = highest temperature [K or Celsius]
T_l = lowest temperature [SAME UNIT AS T_h!!!]
Area = area over which radiation is emitted [m^2]
t = time [s]
emissivity = emissivity of body [-] 0<emisivity<1
absorptivity = absorptivity of body/material [-] 0<absorptivity<1

Output:
Q_rad = radiated heat over time t [J]
"""    

def HeatTransferRadiationEmitter(T_h,T_l,Area,t,emissivity):
        
    Stefan_Boltzmann = 5.670373*10**-8
    Q_radEm = Stefan_Boltzmann*(T_h**4-T_l**4)*Area*t*emissivity
    return Q_radEm
    
def HeatTransferRadiationAbsorber(T_h,T_l,Area,t,absorptivity):
        
    Stefan_Boltzmann = 5.670373*10**-8
    Q_radAb = Stefan_Boltzmann*(T_h**4-T_l**4)*Area*t*absorptivity
    return Q_radAb

"""
Material temperature change

Inputs:
T_init = initial temperature [K or Celsius]
specheat = specific heat of material/system [J/kg*K or J/kg*Celsius]
Q = heat flown into or out of material/system in period t [J]
mass = mass of the system/material [kg]

Output:
T_end = temperature of system after time t [Same unit as T_init]

"""

def TemperatureMaterial(T_init,specheat,Q,mass):
    
    T_end = T_init + (Q/(specheat*mass))
    return T_end

if __name__ == "__main__":
    print HeatTransferConduction(1.4,2,20,10,1,0.003)
    print HeatTransferConvection(2000,1,50,20,1)
    print HeatTransferRadiationEmitter(460,120,1,1,1)
    print HeatTransferRadiationAbsorber(460,120,1,1,1)    
    
