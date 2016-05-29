# -*- coding: utf-8 -*-
"""
Created on Thu May 26 14:09:43 2016

@author: Mathijs

Thermal aircraft
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 26 08:41:42 2016

@author: Mathijs

Thermal control lander
"""

import numpy as np
from HeatTransferFunctions import HeatTransferConduction, HeatTransferConvection, HeatTransferRadiationEmitter, HeatTransferRadiationAbsorber
import Atmosphere


atm = Atmosphere.Atmosphere('preliminary')

height = 3800
latitude = 0
solarLongitude = 0
Temp = atm.temperature(height, latitude, solarLongitude, includeUncertainty=True)

# Temperature in Celsius
T_outside = 200
# Initial equipment temperature
T_eq = 70

mass_Equip = 10.
HeatCapacity_Equip = 10.


R_outerOutshell = 1.3
SpecHeat_Outshell = 10.


R_innerInsul = 0.5
R_outerInsul = 1.
rho_Insul = 5.
SpecHeat_Insul = 10.

R_innerInshell = 0.3
SpecHeat_Inshell = 10.

Q_dotEquip = 10.


# Determine thermal resistances of different layers in the fuselage
TR_Outshell = np.log(R_outerOutshell/R_outerInsul)/(2*np.pi*SpecHeat_Outshell)
TR_Insul = np.log(R_outerInsul/R_innerInsul)/(2*np.pi*SpecHeat_Insul)
TR_Inshell = np.log(R_innerInsul/R_innerInshell)/(2*np.pi*SpecHeat_Inshell)

# Sum thermal resistances to get total resistance
TR_Tot = TR_Outshell+TR_Insul+TR_Inshell

while T_eq < 100:
    dTemp = T_outside-T_eq
    
    dt = 10.
    Q_dot = dTemp/TR_Tot + Q_dotEquip 
    Q_change = Q_dot*dt
    
    T_eq = T_eq + Q_change/(mass_Equip*HeatCapacity_Equip)
    
    print T_eq

    

