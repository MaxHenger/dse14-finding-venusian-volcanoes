# -*- coding: utf-8 -*-
"""
Created on Fri May 27 12:02:26 2016

@author: Mathijs

Aircraft Battery function
"""


'''
Inputs:
Capacity = required battery capacity [Wh]
DOD_system = Depth of Discharge of the battery system [0%<DOD<100%] FOR 80% DOD, WE WILL HAVE ABOUT 1000 CYCLES
SF = safety factor [10% --> SF=0.1, 20% --> SF=0.2, -30% --> SF=-0.3, etc.]

Output:
TotalBatteryWeight = weight of the batteries, including a safety factor [kg]
TotalBatteryVolume = volume of the batteries, including a safety factor [m^3]

Number of cycles is expected to be 500-1000 for 80% DOD, as the 1500 stated in the datasheet is for 30 degrees.
'''

def AircraftBatterySizing(Capacity,DOD,SF):

    # Data of battery
    SpecEnergy = 450. #Wh/kg [See datasheet in dropbox folder: Literature->Battery->OXIS_Li-S_Ultra_Light_Cell_v3.02.pdf]
    SpecDensity = 500. # Wh/L [Same source as SpecEnergy]

    # Define a charging efficiency [Assumed 90%, typical range 80-90% https://blackboard.tudelft.nl/bbcswebdav/pid-2048444-dt-content-rid-7194498_2/courses/28451-131403/Reader%201222%20-%20Spacecraft%20Design%20%28total%29_%202013-2014_v2.pdf]
    ChargeEfficiency = 1.0

    # Define the remaining energy available at mission start (assume 85% of capacity remaining when aircraft is deployed)
    AgeingFactor = 0.85

    # End of cycle remaining charge
    EndOfLife = 0.8

    # Determine the weight of all batteries, and take into account safety factor to account for housing etc.
    TotalBatteryWeight = (Capacity/(SpecEnergy*(DOD/100.)*ChargeEfficiency*AgeingFactor*EndOfLife))*(1.+SF)   # kg
    # Determine the volume of all batteries, and take into account safety factor to account for housing etc.
    TotalBatteryVolume = ((Capacity/(SpecDensity*(DOD/100.)*ChargeEfficiency*AgeingFactor*EndOfLife))/1000.)*(1.+SF)  # m^3

    # Return total battery weight and volume
    return TotalBatteryWeight,TotalBatteryVolume

#if __name__ == '__main__':
    #print  AircraftBatterySizing(50000.,80.,0.)
