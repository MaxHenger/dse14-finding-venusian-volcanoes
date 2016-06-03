# -*- coding: utf-8 -*-
"""
Created on Fri Jun 03 10:45:47 2016

@author: Julius
"""



import numpy as np
import matplotlib.pyplot as plt
import Atmosphere
from ThermalResistanceFunctions import *

atm = Atmosphere.Atmosphere()

class Aircraft:
    def __init__(self,diameter,length,Tinit,SF=1.1,AtmosphericDensity=4.4,BatteryPercentage=0.3,BatteryDensity=1470,internalHeatProd=100,SpecificHeatCapacity=600):
        self.dia = diameter
        self.len = length
        self.vol = np.pi*(self.dia/2.)**2*self.len
        self.rho = BatteryPercentage*BatteryDensity + (1-BatteryPercentage)*AtmosphericDensity
        self.massInternal = self.vol*self.rho
        self.HeatCapacity = self.massInternal*SpecificHeatCapacity
        self.Temp_init=Tinit
        self.insulations=[]
        self.SF=1.1
        
    class Insulation:
        def __init__(self,k,thickness,radius=None):
            self.k=k
            self.t=thickness
            self.r=radius
            self.TR = ThermalResistanceCylinder(self.r+self.t,self.r,self.k)
        
    def addInsulationLayer(self,insu):
        if insu.r==None and len(self.insulations)==0:
            raise AssertionError("Unknown Radius")
        elif insu.r==None:
            insu.r = self.insulations[-1].r+self.insulations[-1].t
        self.insulations.append(insu)
        
    def calc_resistivity(self):
        self.ThermalResis = sum([insu.TR for insu in self.insulations ])
        return self.ThermalResis
        
    def simulate(self,FlightPath=None):
        if FlightPath==None:
            self.FlightPath = self.importFlightPath()
        else:
            self.FlightPath = FlightPath
        #FlightPath = Time, Height, latitude, longitude
        Time = self.FlightPath[0]
        Temp_out = atm.temperature(self.FlightPath[1],self.FlightPath[2],self.FlightPath[3])[1]
        self.calc_resistivity()
        self.Temp_int = [self.Temp_init]
        for i in range(1,len(Time)):
            dt=Time[i]-Time[i-1]
            dTemp = Temp_out[i] - self.Temp_int[i-1]
            q_dot = dTemp/self.ThermalResis *self.SF
            Q_inc = q_dot*dt
            T_inc = Q_inc/self.HeatCapacity
            self.Temp_int.append(self.Temp_int[-1]+T_inc)
            
    def plot_temp(self):
        plt.plot(self.FlightPath[0],self.Temp_int)
        plt.show()

    def importFlightPath(self,FilePath=".\FlightPath.csv"):
        try:
            with open(FilePath) as FP:
                pass
        except IOError:
            pass
        # Import from csv file
        path=np.array([[0,1,2,3],[30000,30001,30002,30003],[0,0,0,0],[0,0,0,0]])
        self.FlightPath=path        
        return path

      
    
    
if __name__=="__main__":
    diameter = 1.1
    length = 4
    Tinit=273
    SF=1.1
    R_insu = 1.
    T_insu = 0.1
    k_insu = 5.8
    ac = Aircraft(diameter,length,Tinit)
    ac.addInsulationLayer(ac.Insulation(k_insu,T_insu,R_insu))
    ac.simulate()
    ac.plot_temp()