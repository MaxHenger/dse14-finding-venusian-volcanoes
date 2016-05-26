# -*- coding: utf-8 -*-
"""
Created on Tue May 10 20:43:17 2016

This file contains the definition of the gravitational model, including all
the J effects.

In general all functions accept and return values in SI units. An exception is
angles, these are entered in degrees. Also attempt to call the functions with
arrays as much as possible to reduce unneccesary calculations.

The normalization is:
http://www.nag.com/numeric/FL/manual/pdf/S/s22aaf.pdf
mathematical proof behind it is found here:
http://en.citizendium.org/wiki/Associated_Legendre_function/Proofs
an example of how the graph should look like is found here:
https://help.scilab.org/docs/5.5.2/en_US/legendre.html

The general use of this file is as following:
    
    import Gravity
    gravE = Gravity.Gravity("Earth")
    gravV = Gravity.Gravity("Venus",accuracy=50)
    gravV(alt,latitude,longitude) # returns gravity at that location
    gravV.a_long(alt,latitude,longitude) # returns longitudinal acceleration at that location
    gravV.a_lat(alt,latitude,longitude) # returns latitudinal acceleration at that location
    
@author: Julius
"""
import numpy as np
import GravityConstants as grav_const
from scipy.special import lpmv,lpmn
#import sympy.mpmath as mp
import matplotlib.pyplot as plt

class Gravity:
    """Gravity class to get the gravity at a specific point"""
    def __init__(self,planet="Venus",method="complex",accuracy=10):
        if planet=="Venus":
            self.constants=grav_const.GravityVenus(method,accuracy)
        elif planet=="Earth":
            self.constants=grav_const.GravityEarth(method,accuracy)
        else:
            raise ValueError("Planet Unknown")
        self.accuracy=accuracy
        self.S=self.constants.S
        self.C=self.constants.C
        self.R=self.constants.RadiusMean
        self.Mu=self.constants.Mu
    
    def __call__(self,altitude,longitude,latitude):
        return self._tinygrav(altitude,longitude,latitude)
        
    def __updateAccuracy__(self,newAccuracy):
        self.accuracy=newAccuracy
        
    def _legn(self,n,x):
        return ( (2*n+1)/2.)**0.5 *lpmv(0,n,x)
        
    def _legnm(self,n,m,x):
        return (-1)**m*( (2*n+1)/2. * float(np.math.factorial(n-m))/np.math.factorial(n+m) )**0.5 * lpmv(m,n,x)
    
    def _updateLegendre(self,latitude):
        #lpmn is not yet normalized
        self.lp,self.derlp = lpmn(self.accuracy,self.accuracy,latitude)
    
    def _normalization(self,n,m):
        return (-1)**m*( (2*n+1)/2. * float(np.math.factorial(n-m))/np.math.factorial(n+m) )**0.5
        
    def _getLP(self,n,m):
        return self._normalization(n,m)*self.lp[m][n]
    
    def _getDerLP(self,n,m):
        return self._normalization(n,m)*self.derlp[m][n]
    
    def __tinygravOld__(self,altitude,longitude,latitude):
        r = altitude+self.R
        g = self.Mu/(r**2) * (1+sum([ (n+1)*(self.R/r)**n * sum([ self._legnm(n,m,np.sin(np.deg2rad(latitude)))*(self.C[n][m]*np.cos(m*np.deg2rad(longitude)) + self.S[n][m]*np.sin(m*np.deg2rad(longitude)) )     for m in range(0,n+1)]) for n in range(2,len(self.C))  ]) )
        return g
    
    def a_longOLD(self,altitude,longitude,latitude):
        r = altitude+self.R
        a_long = self.Mu/(r**2*np.sin(np.deg2rad(latitude)))* sum([ sum([ (self.R/r)**n * self._legnm(n,m,np.cos(np.deg2rad(latitude)))*m*(-self.C[n][m]*np.sin(m*np.deg2rad(longitude)) + self.S[n][m]*np.cos(m*np.deg2rad(longitude)) )  for m in range(0,n)]) for n in range(2,len(self.C)) ])
        return a_long
        
    def _tinygrav(self,altitude,longitude,latitude):
        self._updateLegendre(np.sin(np.deg2rad(latitude)))
        r = altitude+self.R
        g = self.Mu/(r**2) * (1+sum([ (n+1)*(self.R/r)**n * sum([ self._getLP(n,m)*(self.C[n][m]*np.cos(m*np.deg2rad(longitude)) + self.S[n][m]*np.sin(m*np.deg2rad(longitude)) )     for m in range(0,n+1)]) for n in range(2,len(self.C))  ]) )
        return g
        
    def a_lat(self,altitude,longitude,latitude):
        self._updateLegendre(np.cos(np.deg2rad(latitude)))
        r = altitude+self.R
        a_lat = self.Mu/(r**2) * sum([  sum([ (self.R/r)**n * self._getDerLP(n,m)*(self.C[n][m]*np.cos(m*np.deg2rad(m*longitude)) + self.S[n][m]*np.sin(m*np.deg2rad(m*longitude)) ) \
        for m in range(0,n)]) for n in range(2,len(self.C)) ])
        return a_lat
    
    def a_long(self,altitude,longitude,latitude):
        self._updateLegendre(np.cos(np.deg2rad(latitude)))
        r = altitude+self.R
        a_long = self.Mu/(r**2*np.sin(np.deg2rad(latitude)))* sum([ sum([ (self.R/r)**n * self._getLP(n,m)*m*(-self.C[n][m]*np.sin(m*np.deg2rad(longitude)) + self.S[n][m]*np.cos(m*np.deg2rad(longitude)) )  for m in range(0,n)]) for n in range(2,len(self.C)) ])
        return a_long
    
def graph_leg(n=5):
    x=np.arange(-1,1,0.01)
    for m in range(0,n+1):
        y0=lambda x: lpmv(m,n,x)
        yMeg=lambda x: np.sqrt(np.math.factorial(n+m)/(2- (1 if m==0 else 0) )*(2*n+1)*np.math.factorial(n-m) )*lpmv(m,n,x)
        ynorm=lambda x: lpmv(m,n,x) if m==0 else ( (2*n+1)/2. * float(np.math.factorial(n-m))/np.math.factorial(n+m) )**0.5 *lpmv(m,n,x)
        ynormM=lambda x: np.sqrt( (2*n+1)/2. * np.math.factorial(n-m)/np.math.factorial(n+m) ) *lpmv(m,n,x)
        yself=lambda x: grav._legnm(n,m,x) 
        plt.plot(x,yself(x),label="m="+str(m))
        plt.grid(True,which=u'major')
        plt.grid(True,which=u'minor',color="r")
        plt.legend()
    
if __name__=="__main__":   
    grav=Gravity(accuracy=20)
    print(grav(0,30,30))
    print(grav.__tinygravOld__(0,30,30))
    print(grav.a_lat(100000,30,30))
    print(grav.a_long(100000,30,30))
    print(grav.a_longOLD(100000,30,30))
