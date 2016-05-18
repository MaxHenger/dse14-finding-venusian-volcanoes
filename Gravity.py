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
    gravE = Gravity.Gravity("Earth")
    gravV = Gravity.Gravity("Venus",accuracy=50)
    gravV(alt,latitude,longitude)
    
@author: Julius
"""
import numpy as np
import GravityConstants as grav_const
from scipy.special import lpmv,lpmn
import random

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
        return self.__tinygrav__(altitude,longitude,latitude)
        
    def __updateAccuracy__(self,newAccuracy):
        self.accuracy=newAccuracy
        
    def __P1__(self,n,x):
        return lpmv(0,n,x)
        #return 1./ ( (-2)**n * np.math.factorial(n) ) * mp.diff( lambda z: (1-z**2)**n ,x,n)
        
    def __P2__(self,n,m,x):
        return ( (2*n+1)/2. * float(np.math.factorial(n-m))/np.math.factorial(n+m) )**0.5 *lpmv(m,n,x)
        # 
        #return self.__P1__(n,x) if m==0 else  (-1)**m/(2**n*np.math.factorial(n))*(1-x**2)**(m/2.)*mp.diff( lambda x: (x**2-1)**n ,x,n+m)
    
    def __tinygrav__(self,altitude,longitude,latitude):
        r = altitude+self.R
        g = self.Mu/(r**2) * (1+sum([ (n+1)*(self.R/r)**n * sum([ self.__P2__(n,m,np.sin(np.deg2rad(latitude)))*(self.C[n][m]*np.cos(m*np.deg2rad(longitude)) + self.S[n][m]*np.sin(m*np.deg2rad(longitude)) )     for m in range(0,n+1)]) for n in range(2,len(self.C))  ]) )
        return g
        
    def a_lat(self,altitude,longitude,latitude):
        r = altitude+self.R
        a_l = self.Mu/(r**2) * sum([  sum([ (self.R/r)**n  for m in range(0,n)]) for n in range(2,len(self.C)) ])
        return a_l


def test__P1__():
    grav=Gravity()
    import scipy.special
    for degree in range(0,150):
            x=random.random()
            if abs(grav.__P1__(degree,x)-scipy.special.lpn(degree,x))<10**-5:
                print("Pass: ", degree)
            else:
                raise ValueError("Did not pass: "+str(degree))
                
                
def test__P2__():
    grav=Gravity()
    import scipy.special
    for degree in range(0,150):
        for order in range(0,degree+1):
            x=random.random()
            if abs(grav.__P2__(degree,order,x)-scipy.special.lpmv(order,degree,x))<10**-5:
                print(grav.__P2__(degree,order,x),scipy.special.lpmv(order,degree,x))
                print("Pass: ", degree,order)
            else:
                print(grav.__P2__(degree,order,x),scipy.special.lpmv(order,degree,x))
                raise ValueError("Did not pass: "+str(degree)+" "+str(order) )
        
grav=Gravity(accuracy=20)
print(grav(0,0,0))
    