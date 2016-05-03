# -*- coding: utf-8 -*-
"""
Created on Tue May 03 09:55:16 2016

@author: Julius
"""
R_gas =  8.3144598
class gas():
    def __init__(self,volume0, pressure0, temperature0, molar_mass,gamma):
        self.V0=volume0
        self.p0=pressure0
        self.T0=temperature0
        
        self.volume=self.v0
        self.pressure=self.p0
        self.temperature=self.T0
        self.gamma=gamma
        self.molar_mass=molar_mass
        self.Rsp=R_gas/self.molar_mass
        self.mass=self.pressure*self.volume/self.temperature/self.Rsp
        self.moles=self.mass/self.molar_mass
        
    def isentropic(self,newT=None,newP=None,newV=None):
        if type(newT)!=type(None):
            self.volume = self.volume*(newT/self.temperature)**(-1/(self.gamma-1))
            self.pressure = self.pressure*(newT/self.temperature)**(self.gamma/(self.gamma-1))
            self.temperature = newT
        if type(newP)!=type(None):
            self.volume = self.volume*(newP/self.pressure)**(-1/(self.gamma))
            self.temperature = self.temperature*(newP/self.pressure)**((self.gamma-1)/self.gamma)
            self.pressure = newP
        if type(newV)!=type(None):
            self.temperature = self.temperature*(newV/self.volume)**(-(self.gamma-1))
            self.pressure = self.pressure*(newV/self.volume)**(-self.gamma)
            self.volume=newV
    
    def isochoric(self,newT=None,newP=None):
        if type(newT)!=type(None):
            self.pressure=self.pressure/self.temperature*newT
            self.temperature=newT
        if type(newP)!=type(None):
            self.temperature=(self.pressure/self.temperature/newP)**-1
            self.pressure=newP      
            
    def isobaric(self,newT=None,newV=None):
        if type(newT)!=type(None):
            self.volume=self.volume/self.temperature*newT
            self.temperature=newT
        if type(newV)!=type(None):
            self.temperature=(self.volume/self.temperature/newV)**-1
            self.volume=newV
            
    def isothermal(self,newP=None,newV=None):
        if type(newP)!=type(None):
            self.volume=self.pressure*self.volume/newP
            self.pressure=newP
            return True
        if type(newV)!=type(None):
            self.pressure=self.pressure*self.volume/newV
            self.volume=newV
            return True
        return False
            
hydrogen=gas()