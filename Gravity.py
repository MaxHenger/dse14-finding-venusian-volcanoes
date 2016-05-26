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
from scipy.special import lpmv,lpmn,lpn
import random
import sympy.mpmath as mp
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
        return self.__tinygrav__(altitude,longitude,latitude)
        
    def __updateAccuracy__(self,newAccuracy):
        self.accuracy=newAccuracy
        
    def __P1__(self,n,x):
        #return lpmv(0,n,x)
        return 1./ ( (-2)**n * np.math.factorial(n) ) * mp.diff( lambda z: (1-z**2)**n ,x,n,dx=0.00001)
        
    def __P2__(self,n,m,x):
        #return self.__P1__(n,x) if m==0 else ( (2*n+1)/2. * float(np.math.factorial(n-m))/np.math.factorial(n+m) )**0.5 *lpmv(m,n,x)
        return self.__P1__(n,x) if m==0 else (  float(np.math.factorial(n-m))/np.math.factorial(n+m) )**0.5 * (-1)**m/(2**n*np.math.factorial(n))*(1-x**2)**(m/2.)*mp.diff( lambda x: (x**2-1)**n ,x,n+m)
    
    def _updateLegendre(self,latitude):
        self.lp,self.derlp = lpmn(self.accuracy,self.accuracy,latitude)
    
    def normalization(self,n,m):
        delta = 1 if m==0 else 0
        return np.sqrt(np.math.factorial(n+m) / ((2-delta)*(2*n+1)*np.math.factorial(n+m) ) )
    
    def _getLP(self,n,m):
        return self.normalization(n,m)*self.lp[n][m]
    
    def __tinygravOLD__(self,altitude,longitude,latitude):
        self._updateLegendre(latitude)
        r = altitude+self.R
        g = self.Mu/(r**2) * (1+sum([ (n+1)*(self.R/r)**n * sum([ self.__P2__(n,m,np.sin(np.deg2rad(latitude)))*(self.C[n][m]*np.cos(m*np.deg2rad(longitude)) + self.S[n][m]*np.sin(m*np.deg2rad(longitude)) )     for m in range(0,n+1)]) for n in range(2,len(self.C))  ]) )
        return g
        
    def __tinygrav__(self,altitude,longitude,latitude):
        self._updateLegendre(latitude)
        r = altitude+self.R
        g = self.Mu/(r**2) * (1+sum([ (n+1)*(self.R/r)**n * sum([ self.lp[m][n]*self.normalization(n,m)*(self.C[n][m]*np.cos(m*np.deg2rad(longitude)) + self.S[n][m]*np.sin(m*np.deg2rad(longitude)) )     for m in range(0,n+1)]) for n in range(2,len(self.C))  ]) )
        return g
        
    def a_lat(self,altitude,longitude,latitude):
        self._updateLegendre(np.cos(np.deg2rad(latitude)))
        r = altitude+self.R
        a_lat = self.Mu/(r**2) * sum([  sum([ (self.R/r)**n * self.derlp[m][n]*self.normalization(n,m)*(self.C[n][m]*np.cos(m*np.deg2rad(m*longitude)) + self.S[n][m]*np.sin(m*np.deg2rad(m*longitude)) ) \
        for m in range(0,n)]) for n in range(2,len(self.C)) ])
        return a_lat

    def a_long(self,altitude,longitude,latitude):
        self._updateLegendre(np.cos(np.deg2rad(latitude)))
        r = altitude+self.R
        a_long = self.Mu/(r**2*np.sin(np.deg2rad(latitude)))* sum([ sum([ (self.R/r)**n * self.lp[m][n]*m*self.normalization(n,m)*(-self.C[n][m]*np.sin(m*np.deg2rad(longitude)) + self.S[n][m]*np.cos(m*np.deg2rad(longitude)) )  for m in range(0,n)]) for n in range(2,len(self.C)) ])
        
        return a_long
    
def test__P1__():
    grav=Gravity()
    for degree in range(0,100):
            x=random.random()
            if abs(grav.__P1__(degree,x)-lpn(degree,x))<10**-5:
                print("Pass: ", degree)
            else:
                raise ValueError("Did not pass: "+str(degree))
                
                
def test__P2__():
    grav=Gravity()
    for degree in range(0,150):
        for order in range(0,degree+1):
            x=random.random()
            if abs(grav.__P2__(degree,order,x)-lpmv(order,degree,x))<10**-5:
                print(grav.__P2__(degree,order,x),lpmv(order,degree,x))
                print("Pass: ", degree,order)
            else:
                print(grav.__P2__(degree,order,x),lpmv(order,degree,x))
                raise ValueError("Did not pass: "+str(degree)+" "+str(order) )
    
def compare():
    accuracy=20
    grav=Gravity(accuracy=accuracy)
    x=random.random()
    print(x)
    grav._updateLegendre(x)
    for degree in range(0,accuracy):
        for order in range(0,degree+1):
            print(grav.__P2__(degree,order,x),grav._getLP(order,degree))

def compute(x):
    for n in range(0,10):
        for m in range(0,n+1):
            print("n= "+str(n)+"  m= "+str(m))
            print("Numpy: "+str(np.polynomial.legendre.legval2d(m,n,x)))
            print("Scipy: "+str(lpmv(m,n,x))+"\n")
def graph(n=5):
    for m in xrange(0,n):
        x=np.arange(-1,1,0.01)
        y0=lambda x: lpmv(m,n,x)
        yMeg=lambda x: np.sqrt(np.math.factorial(n+m)/(2- (1 if m==0 else 0) )*(2*n+1)*np.math.factorial(n-m) )*lpmv(m,n,x)
        ynorm=lambda x: lpmv(m,n,x) if m==0 else ( (2*n+1)/2. * float(np.math.factorial(n-m))/np.math.factorial(n+m) )**0.5 *lpmv(m,n,x)
        ynormM=lambda x: np.sqrt( (2*n+1)/2. * np.math.factorial(n-m)/np.math.factorial(n+m) ) *lpmv(m,n,x)
        yself=lambda x: grav.__P2__(n,m,x)        
        plt.plot(x,yself(x),label="m="+str(m))
        plt.grid(True,which=u'major')
        plt.grid(True,which=u'minor',color="r")
        plt.legend()
    
    
grav=Gravity(accuracy=20)
print(grav(0,0,30))
print(grav.__tinygravOLD__(0,0,30))
print(grav.a_lat(100000,0,0))
print(grav.a_long(100000,0,0))

graph(5)
#compute(0.5)
#compare()