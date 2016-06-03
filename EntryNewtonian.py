# -*- coding: utf-8 -*-
"""
Created on Thu Jun 02 16:23:00 2016

@author: Julius
"""

import numpy as np
import matplotlib.pyplot as plt

class Newtonian():
    def __init__(self,points,Mach):
        self.geometry(points)
        self.gam = 1.2941
        self.Mach=Mach
        self.CP_max = self.CP(self.Mach,self.gam)
        
    def CP(self,Mach,gamma):
        return ( (self.pfrac(Mach,gamma)) -1 )/( (gamma/2.)*Mach**2)
    
    def pfrac(self,Mach,gamma):
        return (2*gamma*Mach**2 - (gamma-1))/(gamma+1)
        
    def gamma(self,rho0,rho1,Mach):
        n=rho1/rho0
        self.gam = (n+1)/(n-1) - 2*n/( (n-1)**Mach**2 )
        
    def CA(self,theta,alpha):
        return self.CP_max*np.sin(theta)**2 + np.sin(alpha)**2*(1-3*np.sin(theta)**2)
        
    def CN(self,theta,alpha):
        return np.cos(theta)**2*np.sin(2*alpha)
        
    def CM(self,theta,alpha):
        return -2./3*(self.CN(theta,alpha)/(np.tan(theta)*np.cos(theta)**2 ))
        
    def analyse(self,alpha,mode="deg"):
        if mode=="deg":
            alpha=np.deg2rad(alpha)
#        for i in range(1,len(self.points[0])-1):
#            print i
#            self.theta[i]
#            CAi=self.CA(self.theta[i],alpha)*(self.points[1][i]**2-self.points[1][i-1]**2)/self.points[1][-1]**2
#            print CAi
        
        self.CA_T = sum([ \
            self.CA(self.theta[i],alpha)*(self.points[1][i]**2-self.points[1][i-1]**2)/self.points[1][-1]**2 \
            for i in range(0,len(self.points[0])-1)])
        
        self.CN_T = sum([ \
            self.CN(self.theta[i],alpha)*(self.points[1][i]**2-self.points[1][i-1]**2)/self.points[1][-1]**2 \
            for i in range(0,len(self.points[0])-1)])
        
        self.CM_T = self.CM(self.theta[0],alpha) +  sum([ \
            self.CM(self.theta[i],alpha)*(self.points[1][i]**3-self.points[1][i-1]**3)/self.points[1][-1]**3 \
            -self.CN(self.theta[i],alpha)*(self.points[1][i]**2-self.points[1][i-1]**2)*(self.s[i]-self.points[0][i])/self.points[1][-1]**3 \
            for i in range(1,len(self.points[0])-1)])
                
        return self.CA_T, self.CN_T, self.CM_T
                
                
    def geometry(self,points):
        self.points=np.array(points)
        self.x=self.points[0]
        self.y=self.points[1]
        self.dx=self.x[1:]-self.x[:-1]
        self.dy=self.y[1:]-self.y[:-1]
        self.theta=np.tan( self.dy/self.dx )
        self.s = self.y[1:]/(self.dy/self.dx)
        
    def show(self):
        plt.plot(self.x,self.y)
        plt.plot(self.x,-self.y)
        plt.show()
        
    def CAa_plot(self):
        a = np.arange(-5,10,0.1)
        CA=[]
        for i in a:
            CA.append(self.analyse(i)[0])
        plt.plot(a,CA)
        plt.show()
        
    def CNa_plot(self):
        a = np.arange(-5,10,0.1)
        CN=[]
        for i in a:
            CN.append(self.analyse(i)[1])
        plt.plot(a,CN)
        plt.show()
        
    def CMa_plot(self):
        a = np.arange(-5,10,0.1)
        CM=[]
        for i in a:
            CM.append(self.analyse(i)[2])
        plt.plot(a,CM)
        plt.show()
        
        
def test_shield():
    
    y = lambda x: 3*x
    x = np.arange(0,2+0.1,0.1)
    points=[x,y(x)]
    Mach = 10.
    shield = Newtonian(points,Mach)

    return shield
    
    
if __name__=="__main__":
    test=test_shield()
    test.analyse(10)
    #test.show()
    print test.CA_T
    print test.CN_T
    print test.CM_T
    test.CAa_plot()
    