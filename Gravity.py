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
    grav = Gravity.Gravity(Venus)
    grav(alt,latitude,longitude)
    
@author: Julius
"""


import numpy as np
import sympy
import Utility as util
import GravityConstants as grav_const

class Gravity:
    """Gravity class to get the gravity at a specific point"""
    def __init__(self,planet="Venus"):
        if planet=="Venus":
            print 
            self.constants=grav_const.GravityVenus()
            self.J=self.constants.J
            self.lam=self.constants.lam
            self.R=self.constants.RadiusMean
        else:
            raise ValueError("Planet Unknown")
        
    def __call__(self,altitude,longitude,latitude):
        return self.__gravity__(altitude,longitude,latitude)
    def __P1__(self,n,x):
        z=sympy.Symbol("z")
        y = (1-z**2)**n
        yprime = y.diff(z,n)
        func = sympy.lambdify(z,yprime,'numpy')
        value = 1./ ( (-2)**n * util.factorial(n) ) * func(x)
        return value
        
    def __P2__(self,n,m,x):
        z=sympy.Symbol("z")
        y = self.__P1__(n,z)
        yprime = y.diff(z,m)
        func = sympy.lambdify(z,yprime,'numpy')
        value = (1.-x**2)**(m/2.) * func(x) # * diff m self.__P1__(n,x)
        return value
    
    def __gravity__(self,altitude,longitude,latitude):
        r = altitude+self.R
        g = self.constants.Mu/(r)**2 * ( 1 \
-sum([ self.J[n] * (self.R/r)**n * self.__P1__(n,np.sin(np.deg2rad(latitude))) for n in range(2,len(self.J))]) \
+sum([ sum([ self.J[n][m]*(self.R/r)**n*self.__P2__(n,m,np.sin(np.deg2rad(latitude)))*np.cos(m*(longitude-self.lam[n,m])) \
    for m in range(1,n)]) for n in range(2,len(self.J))  ]) )
        return g
        
    def __testJ__(self):
        r = self.R
        lat=np.arange(0,90,15)
        for latitude in lat:
            print(str(latitude)+": ", [ self.J[n] * (self.R/r)**n * self.__P1__(n,np.sin(np.deg2rad(latitude))) for n in range(2,len(self.J))]   )
        