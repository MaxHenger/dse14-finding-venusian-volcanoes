# -*- coding: utf-8 -*-
"""
Created on Fri May 27 12:02:26 2016

@author: Mathijs

Aircraft Battery function
"""

import numpy as np

'''
Inputs:
Capacity = required battery capacity [Wh]
TempBat = Average temperature at which the battery will operate [Degrees Celsius] NOTE: HAS TO BE BETWEEN 80 AND 125 !!!!!
DOD_system = Depth of Discharge of the battery system [0%<DOD<100%]

Output:
Battery weight [kg] and volume [m^3] 

NOTE THIS CODE IS REASONABLY ACCURATE BUT MAY REQUIRE MORE INPUTS
'''

def AircraftBatterySizing(Capacity,TempBat,DOD):
    
    # Battery characteristics from SAFT VL32600-125 battery [SEE DATASHEET IN DROPBOX Folder: Literature->Battery-->VL32600_125.pdf]
    BatteryEnergy = 16.2    # Wh
    BatteryWeight = 0.139   # kg
    Diameter = 0.03205  # m
    Height = 0.06185    # m
    
    # Set safety factors for weight and volume of battery system to account for housing etc.
    SF_weight = 1.05
    SF_volume = 1.05
    
    # Determine volume, specific energy and density of battery
    Volume = 0.25*np.pi*(Diameter**2)*Height    # m^3   
    SpecEnergy = BatteryEnergy/BatteryWeight    # Wh/kg
    Density = BatteryWeight/Volume              # kg/m^3


    # For 100% DOD at 125 degrees we have cycle life of 30 and for 25% DOD a life of 200 cycles. Assume linear scaling we have per % DOD reduction
    Scale = (200./30.)/75.
    # Determine factor of increase in lifetime 
    LifetimeFactor = Scale*(100-DOD)

    # Set up data from 100% DOD
    TempArray = [80,115,125]
    CycleArray = [300,45,30]
    # Set up interpolant
    Parabola = np.polyfit(TempArray,CycleArray,2)

    # Determine number of cycles until capacity is only 70% of original value for 100% DOD
    Cycles = (Parabola[0]*(TempBat**2)+Parabola[1]*TempBat+Parabola[2])*LifetimeFactor
    
    # Define a charging efficiency [Assumed 95%]
    ChargeEfficiency = 0.95
    # Define the remaining energy available after 200 cycles [value used from battery catalogue]
    CapacityReductionFactor = 0.7

    # Determine the weight of all batteries, and take into account safety factor to account for housing etc.
    TotalBatteryWeight = (Capacity/(SpecEnergy*DOD*ChargeEfficiency*CapacityReductionFactor)) * SF_weight    # kg
    # Determine the volume of all batteries, and take into account safety factor to account for housing etc.
    TotalBatteryVolume = (TotalBatteryWeight/Density) * SF_volume # m^3
    
    # Return total battery weight and volume
    return TotalBatteryWeight,TotalBatteryVolume, Cycles
    

if __name__ == '__main__':
    print  AircraftBatterySizing(2500,100,25)