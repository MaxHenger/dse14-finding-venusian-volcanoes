"""
Created on Thu May 26 08:41:42 2016

@author: Mathijs

Thermal performance lander
"""

import numpy as np
import Atmosphere
import matplotlib.pyplot as plt

# Set up atmosphere for the temperature
atm = Atmosphere.Atmosphere('preliminary')

# Initial equipment temperature [Assume initial temperature is temperature of 70 degrees Celsius]
T_eq = 70.      # Celsius
# Initial temperature outer edge insulation [K] and 
T_1 = 70        # Celsius
T_2 = 70        # Celsius

# Set maximum temperature of scientific equipment
T_maxEq = 110.      # Celsius

# Safety factor for the possible connnections running through the insulation and possible holes inside the two shells required for instruments and radiation transfer
SF = 1.1

# Equipment weight
mass_Equip = 20.            # kg
# Heat capacity of combined scientific equipment
HeatCapacity_Equip = 600.   # J/kg*K

R_outerOutshell = 0.45+0.0006   # m
SpecHeat_Outshell = 5.8         # m

# Molar mass of gas used to pressurise inner core
MolarMass = 40.        # g/mol
R_outerInsul = 0.45             # m
R_innerInsul = 0.25+0.000661    # m

# Inner PMC shell characteristics
R_innerInshell = 0.25       # m
SpecHeat_Inshell = 0.381    # m

# Heat generated by the scientific equipment
Q_dotEquip = 20.    # W

# Determine thermal resistances of different layers in the fuselage
TR_Outshell = (R_outerOutshell-R_outerInsul)/(4*np.pi*SpecHeat_Outshell*R_outerOutshell*R_outerInsul)
TR_Inshell = (R_innerInsul-R_innerInshell)/(4*np.pi*SpecHeat_Inshell*R_innerInsul*R_innerInshell)

# Create Empty lists to store the temperatures and time
Tlist = []
T1list = []
T2list = []
timelist = []

# Stepsize
dt = 0.1    # s

# Set up arrays for the system that could be replaced with actual descent data [MAKE SURE DATA HAS TIMESTEP OF 0.1 SEC!!!! FOR ACCURACY]
TimeDescendArray = np.arange(0,3600+dt,dt)
AltitudeArray = np.linspace(38000,0,len(TimeDescendArray))
LatitudeArray = np.zeros(len(TimeDescendArray))
SolarLongitudeArray = np.zeros(len(TimeDescendArray))

for i in xrange(len(TimeDescendArray)):
    
    # Determine temperature at altitude, and position [K]
    Temp = atm.temperature(AltitudeArray[i],LatitudeArray[i],SolarLongitudeArray[i], includeUncertainty=True)

    # Temperature in Celsius on the outside
    T_outside = Temp[2]-273.15      # Celsius

    # Determine average temperature in insulation material 
    T_avgInsul = (T_1 + T_2)/2.  # Celsius    
    SpecHeat_Insul = 0.016*((MolarMass/40.)**-0.35) + 0.000016*T_avgInsul       # W/K*m
    TR_Insul = (R_outerInsul-R_innerInsul)/(4*np.pi*SpecHeat_Insul*R_outerInsul*R_innerInsul)
    
    # Sum thermal resistances to get total resistance
    TR_Tot = TR_Outshell+TR_Insul+TR_Inshell        # Celsius
    
    # Determine temperature differential
    dTemp = T_outside-T_eq      # Celsius
    
    # Determine thermal flow per second
    Q_dot = (dTemp/TR_Tot)*SF + Q_dotEquip   # Joule/s
    # Determine thermal flow over period dt
    Q_change = Q_dot*dt     # Joule
    
    # Determine temperature of the equipment onboard 
    T_eq = T_eq + Q_change/(mass_Equip*HeatCapacity_Equip)      # Celsius
    
    # Determine temperature on inner edge and outer edge of insulation
    dT1 = Q_dot*TR_Outshell     # Celsius
    dT2 = Q_dot*TR_Inshell      # Celsius
    
    # Determine temperature on outer (T_1) and inner edge (T_2) of the insulation
    T_1 = T_1 + dT1     # Celsius
    T_2 = T_eq + dT2    # Celsius
    
    # Fill in temperatures to create plots
    Tlist.append(T_eq)
    T1list.append(T_1)
    T2list.append(T_2)
    
    # Check that equipment has not yet overheated before reaching the surface
    if T_eq > T_maxEq:
        print "Max temperature reached at t", TimeDescendArray[i]
        break


# Reached surface t = time of touchdown
t = TimeDescendArray[-1]

# Make arrays into lists for easy appending
TimeList = list(TimeDescendArray)
AltList = list(AltitudeArray)

# Determine temperature at surface
Temp = atm.temperature(AltitudeArray[-1],LatitudeArray[-1],SolarLongitudeArray[-1], includeUncertainty=True)
    
# For as long as the temperature of equipment is below maximum continue computations
while T_eq < T_maxEq:
    
    # Increase time with respect to initial touchdown time
    t = t + dt      # s

    # Temperature in Celsius on the outside
    T_outside = Temp[2]-273.15

    # Determine average temperature of insulation
    T_avgInsul = (T_1 + T_2)/2.  # Celsius    
    # Determine insulation performance and Thermal resistance
    SpecHeat_Insul = 0.016*((MolarMass/40.)**-0.35) + 0.000016*T_avgInsul       # W/K*m
    TR_Insul = (R_outerInsul-R_innerInsul)/(4*np.pi*SpecHeat_Insul*R_outerInsul*R_innerInsul)
    
    # Sum thermal resistances to get total resistance
    TR_Tot = TR_Outshell+TR_Insul+TR_Inshell        # Celsius
    
    # Determine temperature differential
    dTemp = T_outside-T_eq      # Celsius
    
    # Determine thermal flow per second
    Q_dot = dTemp/TR_Tot + Q_dotEquip   # Joule/s
    # Determine thermal flow over period dt
    Q_change = Q_dot*dt     # Joule
    
    # Determine temperature of the equipment onboard 
    T_eq = T_eq + Q_change/(mass_Equip*HeatCapacity_Equip)      # Celsius
    
    # Determine temperature on inner edge and outer edge of insulation
    dT1 = Q_dot*TR_Outshell     # Celsius
    dT2 = Q_dot*TR_Inshell      # Celsius
    
    # Determine temperature on outer (T_1) and inner edge (T_2) of the insulation
    T_1 = T_1 + dT1     # Celsius
    T_2 = T_eq + dT2    # Celsius
    
    # Fill in temperatures to create plots
    Tlist.append(T_eq)
    T1list.append(T_1)
    T2list.append(T_2)
    TimeList.append(t)
    AltList.append(0)

# Create plots
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(TimeList, AltList,'g-')
ax2.plot(TimeList, Tlist,'b-')
ax2.plot(TimeList, T1list, 'r-')

ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Altitude (m)')
ax2.set_ylabel('Temperature (Celsius)')
    
plt.grid()
plt.show()

print "Time of Mission End (s)"
print TimeList[-1]
print "Time spent on surface (s)"
print TimeList[-1]-TimeDescendArray[-1]