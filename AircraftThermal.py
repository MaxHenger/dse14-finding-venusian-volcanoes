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

class AircraftThermal:
    def __init__(self,diameter,length,SF=1.1,AtmosphericDensity=4.4,BatteryPercentage=0.3,BatteryDensity=1470,internalHeatProd=100,SpecificHeatCapacity=600):
        self.dia = diameter
        self.len = length
        self.vol = np.pi*(self.dia/2.)**2*self.len
        self.rho = BatteryPercentage*BatteryDensity + (1-BatteryPercentage)*AtmosphericDensity
        self.massInternal = self.vol*self.rho
        self.HeatCapacity = self.massInternal*SpecificHeatCapacity
        self.internalPower = internalHeatProd
        self.SF=SF
        self.insulations=[]
        
    class Insulation:
        def __init__(self,k,thickness,radius=None,length=None):
            self.k=k
            self.t=thickness
            self.r=radius
            self.l=length
        def calc_res(self):
            self.TR_c = ThermalResistanceCylinder(self.r+self.t,self.r,self.k,self.l)
            self.TC_s = ThermalResistanceSphere(self.r+self.t,self.r,self.k)
            self.TR = (1/self.TR_c + 1/self.TC_s)**-1
            
    def addInsulationLayer(self,insu):
        if insu.r==None and len(self.insulations)==0:
            raise AssertionError("Unknown Radius")
        elif insu.r==None:
            insu.r = self.insulations[-1].r+self.insulations[-1].t
        if insu.l==None:
            insu.l=self.len
        insu.calc_res()
        self.insulations.append(insu)
        
    def calc_resistivity(self):
        self.ThermalResis = sum([insu.TR for insu in self.insulations ])
        return self.ThermalResis
        
    def simulate(self,FlightPath=None,repetition=1):
        self.repetition=repetition
        if FlightPath==None:
            self.FlightPath = self.importFlightPath()
        else:
            self.FlightPath = FlightPath
        #FlightPath = Time, Height, latitude, longitude
        Time = self.FlightPath[0]
        self.Temp_out = [atm.temperature(self.FlightPath[1][0],0,0)[1]]
        self.Temp_init=self.Temp_out[-1]
        self.Temp_int = [self.Temp_init]
        self.Time_sim = [self.FlightPath[0][0]]
        self.H_sim = [self.FlightPath[1][0]]
        self.calc_resistivity()
        for rep in range(0,repetition):
            for i in range(1,len(Time)):
                dt=Time[i]-Time[i-1]
                self.Time_sim.append(self.Time_sim[-1]+dt)
                self.H_sim.append(self.FlightPath[1][i])
                self.Temp_out.append(atm.temperature(self.H_sim[-1],0,0)[1])
                dTemp = self.Temp_out[i] - self.Temp_int[-1]
                q_dot = dTemp/self.ThermalResis *self.SF + self.internalPower
                Q_inc = q_dot*dt
                T_inc = Q_inc/self.HeatCapacity
                self.Temp_int.append(self.Temp_int[-1]+T_inc)
                
            uptime=60*60
            dt=0.5
            for i in np.arange(0,uptime,dt):
                self.Time_sim.append(self.Time_sim[-1]+dt)
                self.H_sim.append(self.FlightPath[1][-1])
                self.Temp_out.append(atm.temperature(self.H_sim[-1],0,0)[1])
                dTemp = self.Temp_out[-1] - self.Temp_int[-1]
                q_dot = dTemp/self.ThermalResis *self.SF + self.internalPower
                Q_inc = q_dot*dt
                T_inc = Q_inc/self.HeatCapacity
                self.Temp_int.append(self.Temp_int[-1]+T_inc)
                

            
    def plot_sim(self):
        plt.plot(self.Time_sim,np.array(self.Temp_int)-273.15,label="Aircraft")
        plt.plot(self.Time_sim,np.array(self.Temp_out)-273.15,label="Atmosphere")
        plt.grid(True)
        plt.legend(loc=1)
        plt.xlabel(r"Time [s]",fontsize=14)
        plt.ylabel(r"Temperature [$\degree$ C]",fontsize=14)
        plt.tight_layout()
        plt.show()
    
    def plot_temp(self):
        plt.plot(self.Time_sim,np.array(self.Temp_out)-273.15)
        plt.show()
        
    def plot_FP(self):
        plt.plot(self.Time_sim,self.H_sim)
        plt.show()        

    def importFlightPath(self,FilePath=".\data\FlightPath.csv",mode="csv"):
        path=[[],[]]        
        if mode=="csv":
            path = np.genfromtxt(FilePath,delimiter=";")
            path = path.T
        elif mode=="dat":
            import TrackStorage as TS
            stor = TS.DataStorage()
            stor.load(FilePath)
            t=stor.getVariable("timeTotal").getValues()
            h=stor.getVariable("heightTotal").getValues()
            path=[t,h]
        else:
            raise AttributeError("Unknown Mode")
            
        self.FlightPath=path        
        return path

      
    
    
if __name__=="__main__":
    diameter = 1
    length = 3
    SF=1.2
    
    R_insu = 1.
    T_insu = 5./1000
    k_insu = 20*10**-3 # 

    T_comp = 2./1000
    k_comp = 0.530    
    
    T_tef = 0.05/1000
    k_tef = 0.238
    
    internalHeatProd=0#5000
    ac = AircraftThermal(diameter,length,SF=SF,internalHeatProd=internalHeatProd)
    ac.addInsulationLayer(ac.Insulation(k_insu,T_insu,R_insu))
    ac.addInsulationLayer(ac.Insulation(k_comp,T_comp))    
    ac.addInsulationLayer(ac.Insulation(k_tef,T_tef))
    FP = ac.importFlightPath(".\data\stitched_62000.0at3.5to38000.0at38000.0.dat","dat")
    ac.simulate(FP,repetition=3)
    ac.plot_sim()
    #ac.plot_FP()
