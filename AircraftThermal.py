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
        self.Temp_init=self.Temp_out[0]
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
                
    def track_descend(self,trackdown,T_init=None,dis_eff=0.95):
        FP = trackdown
        #FlightPath = Time, Height, Power
        Time_track = FP[0]
        Temp_out = [atm.temperature(FP[1][0],0,0)[1]]
        if T_init==None:        
            Temp_init= Temp_out[0]
        else:
            Temp_init=T_init
        Temp_int = [Temp_init]
        Time_sim = [FP[0][0]]
        H_sim    = [FP[1][0]]
        self.calc_resistivity()
        try:
            internal_power = FP[2][0]
        except TypeError:
            internal_power = FP[2]
            
        for i in range(1,len(Time_track)):
            dt=Time_track[i]-Time_track[i-1]
            Time_sim.append(Time_sim[-1]+dt)
            H_sim.append(FP[1][i])
            Temp_out.append(atm.temperature(H_sim[-1],0,0)[1])
            dTemp = Temp_out[-1] - Temp_int[-1]
            q_dot = dTemp/self.ThermalResis *self.SF + internal_power*(1-dis_eff)
            Q_inc = q_dot*dt
            T_inc = Q_inc/self.HeatCapacity
            Temp_int.append(Temp_int[-1]+T_inc)
        return Time_sim,H_sim,Temp_int
    
    def track_climb(self,trackup,T_init=None,dis_eff=0.95,Accuracy=1,max_T=80+273):
        FP = trackup
        #FlightPath = Time, Height, Power
        Time_track = FP[0]
        Temp_out = [atm.temperature(FP[1][0],0,0)[1]]
        if T_init==None:        
            Temp_init= Temp_out[0]
        else:
            Temp_init=T_init
        Temp_int = [Temp_init]
        Time_sim = [FP[0][0]]
        H_sim    = [FP[1][0]]
        self.calc_resistivity()
        for i in range(1,len(Time_track)):
            dt=Time_track[i]-Time_track[i-1]
            Time_sim.append(Time_sim[-1]+dt)
            H_sim.append(FP[1][i])
            Temp_out.append(atm.temperature(H_sim[-1],0,0)[1])
            dTemp = Temp_out[-1] - Temp_int[-1]
            q_dot = dTemp/self.ThermalResis *self.SF + FP[2][i]*(1-dis_eff)
            Q_inc = q_dot*dt
            T_inc = Q_inc/self.HeatCapacity
            Temp_int.append(Temp_int[-1]+T_inc)
            if Temp_int[-1]>(max_T):
                break
        return Time_sim,H_sim,Temp_int 
        
    def time_down(self,FPdown, FPup, t_down,step, max_T=80+273, T_init=None, dis_eff=0.95,dt=0.5,Accuracy=1):
        #print ("\n") 
        #print ("Step: ",step)
        #print ("Temp Init: ",T_init-273)
        Time_sim,H_sim,Temp_int=self.track_descend(FPdown,T_init)
        #print ("Temp dive: ",Temp_int[-1]-273)
        #print ("Alt down: ",H_sim[-1])
        #print ("Time down: ",t_down)
        downtime=t_down
        dt=dt
        for i in np.arange(0,downtime,dt):
            Time_sim.append(Time_sim[-1]+dt)
            H_sim.append(FPdown[1][-1])
            dTemp = atm.temperature(FPdown[1][-1],0,0)[1] - Temp_int[-1]
            q_dot = dTemp/self.ThermalResis *self.SF + FPdown[2]*(1-dis_eff)
            Q_inc = q_dot*dt
            T_inc = Q_inc/self.HeatCapacity
            Temp_int.append(Temp_int[-1]+T_inc)
            if Temp_int[-1]>=(max_T):
                break
        #print ("Temp down",Temp_int[-1]-273)
        if Temp_int[-1]<max_T:
            Time_sim,H_sim,Temp_int = self.track_climb(FPup,T_init=Temp_int[-1],Accuracy=Accuracy,max_T=max_T)
            print ("Temp climb",Temp_int[-1]-273)
        
        #print ("Delta T: ",Temp_int[-1]-max_T)
        
        if step<Accuracy:
            return t_down
            
        if Temp_int[-1]<max_T:
            return self.time_down(FPdown, FPup, t_down+step/2.,step/2., max_T, T_init, dis_eff,dt)
        elif Temp_int[-1]>max_T:
            return self.time_down(FPdown, FPup, t_down-step/2.,step/2., max_T, T_init, dis_eff,dt,Accuracy)
        
        
        
    
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
            t=stor.getVariable("time").getValues()
            h=stor.getVariable("height").getValues()
            p=stor.getVariable("power").getValues()
            path=[t,h,p]
        else:
            raise AttributeError("Unknown Mode")
            
        #self.FlightPath=path        
        return path
    
    def loadFlightPath(self,Storage):
        stor = Storage
        t=stor.getVariable("time").getValues()
        h=stor.getVariable("height").getValues()
        p=stor.getVariable("power").getValues()
        path=[t,h,p]
        return path
        

def easy_setup():
    diameter = 1
    length = 3
    SF=1.2
    
    R_insu = diameter/2.
    T_insu = 5./1000
    k_insu = 20*10**-3 # 

    T_comp = 2./1000
    k_comp = 0.530    
    
    T_tef = 0.05/1000
    k_tef = 0.238
    
    V = 2*np.pi* R_insu*length*T_insu + np.pi*R_insu**2*T_insu
    rho = 130
    print("Mass: ",V*rho)
    Max_W = 41391
    Dis_eff = 0.95
    internalHeatProd=Max_W*(1-Dis_eff)
    print("Internal Heating: ",internalHeatProd)
    ac = AircraftThermal(diameter,length,SF=SF,internalHeatProd=internalHeatProd)
    ac.addInsulationLayer(ac.Insulation(k_insu,T_insu,R_insu))
    ac.addInsulationLayer(ac.Insulation(k_comp,T_comp))    
    ac.addInsulationLayer(ac.Insulation(k_tef,T_tef))
    #FP = ac.importFlightPath(".\data\stitched_62000.0at3.5to38000.0at38000.0.dat","dat")
    #ac.simulate(FP,repetition=3)
    return ac
    
if __name__=="__main__":
    ac=easy_setup()
    
    t_down = 10000
    step = t_down
    
    ### FP = [Time, Altitude, Power]
    FPdown = ac.importFlightPath(".\data\dive_0.0_66658.299at2.634_32000.0at-21.536.dat","dat")
    FPup = ac.importFlightPath(".\data\climb_0.0_32000.0_to_66658.299.dat","dat")
    time = ac.time_down(FPdown, FPup,t_down,step,T_init=10+273,dt=1,Accuracy=10) 

    #Comments to Max
    #
    # call easy_setup() for the current setuo
    # define your two flight track (use ac.loadFlightPath(storage) )
    # call ac.time_down
    # Accuracy will determine the smallest time step size for the bottom time
    # dt is for the down time time step for the simulation, 1 should suffice
    # T_init is 10 deg, or smthing around that
    
    
    
    
