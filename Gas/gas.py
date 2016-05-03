# -*- coding: utf-8 -*-
"""
Created on Tue May 03 09:55:16 2016

@author: Julius
"""
R_gas =  8.3144598
class gas():
    def __init__(self,volume0, pressure0, temperature0, molar_mass,gamma):
        self.v0=volume0
        self.p0=pressure0
        self.T0=temperature0
        
        self.volume=self.v0
        self.pressure=self.p0
        self.temperature=self.T0
        self.gamma=gamma
        self.molar_mass=molar_mass
        if molar_mass>1:
            self.molar_mass/=1000.
        self.Rsp=R_gas/self.molar_mass
        self.mass=self.pressure*self.volume/self.temperature/self.Rsp
        self.moles=self.mass/self.molar_mass
        self.density= self.pressure/(self.Rsp*self.temperature)

    def ideal(self,newT=None,newP=None,newV=None):
        if type(newT)!=type(None) and type(newP)!=type(None):
            self.volume=self.pressure*self.volume/self.temperature * newT/(newP)
        if type(newT)!=type(None) and type(newV)!=type(None):
            self.pressure=self.pressure*self.volume/self.temperature * newT/(newV)
        if type(newP)!=type(None) and type(newV)!=type(None):
            self.temperature=self.pressure*self.volume/self.temperature /(newP*newV)
            
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
        self.density= self.pressure/(self.Rsp*self.temperature)
    
    def isochoric(self,newT=None,newP=None):
        if type(newT)!=type(None):
            self.pressure=self.pressure/self.temperature*newT
            self.temperature=newT
        if type(newP)!=type(None):
            self.temperature=(self.pressure/self.temperature/newP)**-1
            self.pressure=newP
        self.density= self.pressure/(self.Rsp*self.temperature)
            
    def isobaric(self,newT=None,newV=None):
        if type(newT)!=type(None):
            self.volume=self.volume/self.temperature*newT
            self.temperature=newT
        if type(newV)!=type(None):
            self.temperature=(self.volume/self.temperature/newV)**-1
            self.volume=newV
        self.density= self.pressure/(self.Rsp*self.temperature)
            
    def isothermal(self,newP=None,newV=None):
        if type(newP)!=type(None):
            self.volume=self.pressure*self.volume/newP
            self.pressure=newP
        if type(newV)!=type(None):
            self.pressure=self.pressure*self.volume/newV
            self.volume=newV
        self.density= self.pressure/(self.Rsp*self.temperature)
            
if __name__=="__main__":
    hydrogen=gas(112,0.0369*10**5,229.8,2.016,1.41)
    
    helium=gas(112,0.0369*10**5,229.8,8,1.66)