# -*- coding: utf-8 -*-
"""
Created on Wed Jun 01 23:06:03 2016

@author: Julius
"""
from scipy.interpolate import griddata

class Vehicle:
    def __init__(self,name="Unknown"):
        self.name=name
        
    def geometry(self,m,S,b,c,Ixx,Iyy,Izz,Ixz):
        self.m=m
        self.S=S
        self.b=b
        self.c=c
        self.Ixx=Ixx
        self.Iyy=Iyy
        self.Izz=Izz
        self.Ixz=Ixz
        
    def update(self,alpha,beta,mach):
        self.CDa(alpha,mach)
        self.CDM(alpha,mach)
        self.CLa(alpha,mach)
        self.CLM(alpha,mach)
        self.Cma(alpha,mach)
        self.CmM(alpha,mach)
        
        self.Clb(beta)
        self.Cnb(beta)
        self.CSb(beta)
        
        self.Clda(alpha,mach)
        self.Cmda(alpha,mach)
        self.Cnda(alpha,mach)
        self.Cndr(alpha,mach)
    
    def _grids(self,CD,CL,Cm,Cl,Cn,CS):
        self.CDgrd=CD
        self.CLgrd=CL
        self.Cmgrd=Cm
        self.Clgrd=Cl
        self.Cngrd=Cn
        self.CSgrd=CS
    
    def _der_grd(self,grd,x,y=0,D=1,h=0.001):
        grd
    
    def CDa(self,alpha,mach):
        pass
    def CDM(self,alpha,mach):
        pass
    
    def CLa(self,alpha,mach):
        pass
    def CLM(self,alpha,mach):
        pass

    def Cma(self,alpha,mach):
        pass
    def CmM(self,alpha,mach):
        pass

    def Clb(self,beta):
        pass
    def Cnb(self,beta):
        pass
    def CSb(self,beta):
        pass

    def Clda(self,alpha,mach):
        pass
    def Cmda(self,alpha,mach):
        pass
    def Cnda(self,alpha,mach):
        pass
    def Cndr(self,alpha,mach):
        pass


def Test_vehicle():
    Test_Vehicle=Vehicle("Test")
    ### Geometric stuff
    m=26029. # kg mass
    S=110. # m2 surface area
    b=13.
    c=23
        ### Inertia stuff
    Ixx=1000
    Iyy=1000
    Izz=1000
    Ixz=100

    Test_Vehicle.geometry(m,S,b,c,Ixx,Iyy,Izz,Ixz)
    
    
    
    
    return Test_Vehicle
