# -*- coding: utf-8 -*-
"""
Created on Fri Jun 03 10:45:47 2016

@author: Julius
"""



import numpy as np
import matplotlib.pyplot as plt
import Atmosphere

atm = Atmosphere.Atmosphere()

class Aircraft:
    def __init__(self,diameter,length,Tinit,AtmosphericDensity=4.4,BatteryPercentage=0.3,BatteryDensity=1470,internalHeatProd=100,SpecificHeatCapacity=600):
        self.dia = diameter
        self.len = length
        self.vol = np.pi*(self.dia/2.)**2*self.len
        self.rho = BatteryPercentage*BatteryDensity + (1-BatteryPercentage)*AtmosphericDensity
        self.massInternal = self.vol*self.rho
        self.HeatCapacity = self.massInternal*SpecificHeatCapacity
        self.Temp_init=Tinit
        
    class Insulation:
        def __init__(self,R_outer,R_inner,SpecificHeat,ThermalResist):
            self.R_o = R_outer
            self.R_i = R_inner
            self.K_sp = SpecificHeat
    
    
    def addInsulationLayer(self,thickness,radius,):
        pass


def importFlightPath(FilePath=".\JuliusAwesomeFile.csv"):
    try:
        with open(FilePath) as FP:
            pass
    except IOError:
        pass
    # Import from csv file
    path=np.array([[0,1,2,3],[30000,30001,30002,30003],[0,0,0,0],[0,0,0,0]])
    return path

      
def simulateThermal(Aircraft,FlightPath,SF):
    Time = FlightPath[0]
    Temperature = atm.temperature(FlightPath[1],FlightPath[2],FlightPath[3])[1]
    print Temperature
    Temp_internal = [Aircraft.Temp_init]
    for i in range(1,len(Time)):
        dt=Time[i+1]-Time[i]
        dTemp = Temperature[i] - Temp_internal[i-1]
        
        q_dot = dTemp/Aircraft.insuResist *SF
        Q_inc = q_dot*dt
        T_inc = Q_inc/Aircraft.HeatCapacity
        Temp_internal.append(Temp_internal[-1]+T_inc)
    
    
if __name__=="__main__":
    
    diameter = 1.1
    length = 4
    Tinit=273
    SF=1.1
    ac = Aircraft(diameter,length,Tinit)
    fp = importFlightPath()
    simulateThermal(ac,fp,SF)