# -*- coding: utf-8 -*-
"""
Created on Tue May 10 20:43:17 2016

This file contains the definition of the gravitational model, including all
the J effects.

In general all functions accept and return values in SI units. An exception is
angles, these are entered in degrees. Also attempt to call the functions with
arrays as much as possible to reduce unneccesary calculations.

The general use of this file is as following:
    
    import Gravity
    grav = Gravity.Gravity()
    grav(alt,latitude,longitude)
    
@author: Julius
"""


import numpy as np
import sympy
import Utility as util
import GravityConstants as grav_const
import sympy.mpmath as mp

class Gravity:
    """Gravity class to get the gravity at a specific point"""
    def __init__(self,planet="Venus",GravModel="Magellen",accuracy=10):
        
        if planet=="Venus":
            self.constants=grav_const.GravityVenus(GravModel,accuracy)
            self.GravModel=GravModel
            self.accuracy=accuracy
            if GravModel=="Magellen":
                self.S=self.constants.S
                self.C=self.constants.C
            elif GravModel=="simple":
                self.J=self.constants.J
                self.lam=self.constants.lam
            self.R=self.constants.RadiusMean
        
        elif planet=="Earth":
            self.constants=grav_const.GravityEarth()
            self.J=self.constants.J
            self.lam=self.constants.lam
            self.GravModel="simple"
            self.R=self.constants.RadiusMean
            
        else:
            raise ValueError("Planet Unknown")
        
    def __call__(self,altitude,longitude,latitude):
        return self.__gravity__(altitude,longitude,latitude)
        
    def __P1__(self,n,x):
        return 1./ ( (-2)**n * util.factorial(n) ) * mp.diff( lambda z: (1-z**2)**n ,x,n)
        
    def __P2__(self,n,m,x):
        return (1.-x**2)**(m/2.) * mp.diff( lambda z: 1./ ( (-2)**n * util.factorial(n) ) * (1-z**2)**n ,x,m)
    
    def __gravity__(self,altitude,longitude,latitude):
        r = altitude+self.R
        if self.GravModel=="simple":
            g = self.constants.Mu/(r**2) * ( 1 \
    -sum([ self.J[n][0] * (n+1)*(self.R/r)**n * self.__P1__(n,np.sin(np.deg2rad(latitude))) for n in range(2,len(self.J))]) \
    +sum([ sum([ self.J[n][m]*(n+1)*(self.R/r)**n*self.__P2__(n,m,np.sin(np.deg2rad(latitude)))*np.cos(m*(np.deg2rad(longitude-self.lam[n,m]))) \
        for m in range(1,n)]) for n in range(2,len(self.J))  ]) )
            
        elif self.GravModel=="Magellen":
#            r0 = altitude+self.R
#            r = sympy.Symbol("r")
#            U = -self.constants.Mu/(r) * ( 1 \
#        -sum([ -self.C[n][0] * (self.R/r)**n * self.__P1__(n,np.sin(np.deg2rad(latitude))) for n in range(2,len(self.C))]) \
#        +sum([ sum([ (self.R/r)**n*self.__P2__(n,m,np.sin(np.deg2rad(latitude)))*(self.C[n][m]*np.cos(m*np.deg2rad(longitude)) + self.S[n][m]*np.sin(m*np.deg2rad(longitude)) ) \
#        for m in range(1,n+1)]) for n in range(2,len(self.C))  ]) )
#            Uprime = U[1].diff(r)
#            func = sympy.lambdify(r,Uprime,'numpy')
#            g1 = func(r0)
            
            g = self.constants.Mu/(r**2) * ( 1 \
    -sum([ -self.C[n][0] *(n+1)*(self.R/r)**n * self.__P1__(n,np.sin(np.deg2rad(latitude))) for n in range(2,len(self.C))]) \
    +sum([ sum([  (n+1)*(self.R/r)**n*self.__P2__(n,m,np.sin(np.deg2rad(latitude)))*(self.C[n][m]*np.cos(m*np.deg2rad(longitude)) + self.S[n][m]*np.sin(m*np.deg2rad(longitude)) ) \
    for m in range(1,n+1)]) for n in range(2,len(self.C))  ]) )
        
        else:
            raise ValueError("unknown GravModel")
            
        return g
        
        
    def __testJ__(self):
        r = self.R
        lat=np.arange(0,90,15)
        for latitude in lat:
            print(str(latitude)+": ", [ self.J[n] * (self.R/r)**n * self.__P1__(n,np.sin(np.deg2rad(latitude))) for n in range(2,len(self.J))]   )

if __name__=="__main__":
    grav=Gravity(accuracy=150)
    print(grav(0,0,0))[1]