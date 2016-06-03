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
        #return 2*np.sin(theta)**2 + np.sin(alpha)**2*(1-3*np.sin(theta)**2)
        return self.CP_max*np.sin(theta)**2 + np.sin(alpha)**2*(1-3*np.sin(theta)**2)
        
    def CN(self,theta,alpha):
        #return np.cos(theta)**2*np.sin(2*alpha)
        return self.CP_max/2.*np.cos(theta)**2*np.sin(2*alpha)
        
    def CM(self,theta,alpha):
        #return -2./3*(self.CN(theta,alpha)/(np.tan(theta)*np.cos(theta)**2 ))
        return -self.CP_max/3.*(self.CN(theta,alpha)/(np.tan(theta)*np.cos(theta)**2 ))
        
    def analyse(self,alpha,mode="deg"):
        if mode=="deg":
            alpha=np.deg2rad(alpha)
        
        self.CA_T = sum([ \
            self.CA(self.theta[i],alpha)*(self.points[1][i+1]**2-self.points[1][i]**2)/self.points[1][-1]**2 \
            for i in range(0,len(self.dx))])
        
        self.CN_T = sum([ \
            self.CN(self.theta[i],alpha)*(self.points[1][i+1]**2-self.points[1][i]**2)/self.points[1][-1]**2 \
            for i in range(0,len(self.dx))])
        
        self.CM_T = self.CM(self.theta[0],alpha) +  sum([ \
            self.CM(self.theta[i],alpha)*(self.points[1][i+1]**3-self.points[1][i]**3)/self.points[1][-1]**3 \
            -self.CN(self.theta[i],alpha)*(self.points[1][i+1]**2-self.points[1][i]**2)*(self.s[i]-self.points[0][i+1])/self.points[1][-1]**3 \
            for i in range(1,len(self.dx))])
                
        return self.CA_T, self.CN_T, self.CM_T
                
                
    def geometry(self,points):
        self.points=np.array(points)
        self.x=self.points[0]
        self.y=self.points[1]
        self.dx=self.x[1:]-self.x[:-1]
        self.dy=self.y[1:]-self.y[:-1]
        self.theta=np.arctan( self.dy/self.dx )
        self.s = self.y[1:]/(self.dy/self.dx)
        
    def show(self):
        plt.plot(self.x,self.y,color="b")
        plt.plot(self.x,-self.y,color="b")
        plt.show()
        
    def show_flight(self,alpha):
        a = np.deg2rad(alpha)
        self.T=np.matrix([[np.cos(a),np.sin(a)],[-np.sin(a),np.cos(a)]])
        xUnew=np.zeros(len(self.x))
        yUnew=np.zeros(len(self.y))
        xLnew=np.zeros(len(self.x))
        yLnew=np.zeros(len(self.y))
        for i in range(len(self.x)):
            xUnew[i]=(self.T*np.matrix([[self.x[i]],[self.y[i]]]))[0]
            yUnew[i]=(self.T*np.matrix([[self.x[i]],[self.y[i]]]))[1]
            xLnew[i]=(self.T*np.matrix([[self.x[i]],[-self.y[i]]]))[0]
            yLnew[i]=(self.T*np.matrix([[self.x[i]],[-self.y[i]]]))[1]
        plt.plot(xUnew,yUnew,color="r")
        plt.plot(xLnew,yLnew,color="r")
        maxX=max(xUnew[-1],xLnew[-1])
        maxY=max(abs(yLnew[-1]),abs(yUnew[-1]))
        maxT=max(maxX,maxY)
        plt.xlim((-maxT,maxT))
        plt.ylim((-maxT,maxT))
        plt.show()
        
    def CAa_plot(self):
        a = np.arange(-5,10,0.1)
        CAa=[]
        for i in a:
            CAa.append(self.analyse(i)[0])
        plt.plot(a,CAa)
        plt.show()
        
    def CNa_plot(self):
        a = np.arange(-5,10,0.1)
        CNa=[]
        for i in a:
            CNa.append(self.analyse(i)[1])
        plt.plot(a,CNa)
        plt.show()
        
    def CMa_plot(self):
        a = np.arange(-5,10,0.1)
        CMa=[]
        for i in a:
            CMa.append(self.analyse(i)[2])
        plt.plot(a,CMa)
        plt.show()
        
    def CAM_plot(self,alpha=5):
        M = np.arange(1,20,0.5)
        
        CA=[]
        for i in M:
            self.CP_max = self.CP(i,self.gam)
            CA.append(self.analyse(alpha)[0])
        self.CP_max = self.CP(self.Mach,self.gam)
        plt.plot(M,CA)
        plt.show()
        
    def CNM_plot(self,alpha=5):
        M = np.arange(1,20,0.5)
        CN=[]
        for i in M:
            self.CP_max = self.CP(i,self.gam)
            CN.append(self.analyse(alpha)[1])
        self.CP_max = self.CP(self.Mach,self.gam)
        plt.plot(M,CN)
        plt.show()
        
    def CMM_plot(self,alpha=5):
        M = np.arange(1,20,0.5)
        CM=[]
        for i in M:
            self.CP_max = self.CP(i,self.gam)
            CM.append(self.analyse(alpha)[2])
        self.CP_max = self.CP(self.Mach,self.gam)
        plt.plot(M,CM)
        plt.show()
        
    def CP_plot(self):
        M = np.arange(1,20,0.5)
        CP=[]
        for i in M:
            CP.append(self.CP(i,self.gam))
        self.CP_max = self.CP(self.Mach,self.gam)
        plt.plot(M,CP)
        plt.show()
        
        
def test_shield():
    depth = 3.
    width = 4.6
    depth_nose = 0.1
    width_nose = 1
    def y(x):
        if x <= depth_nose:
            return width_nose*(x/depth_nose)**0.3
        else:
            return width_nose+width*((x-depth_nose)/depth)**0.8
    dt=0.0001
    x = np.arange(0,depth+dt,dt)
    yout=np.zeros(len(x))
    for i,x_i in enumerate(x):
        yout[i]=y(x_i)
    points=np.array([x,yout])
    #print x
    #print yout
    
    #points=[[0,0.1,0.2,0.3,0.4,0.5,2],[0,0.3162,0.447,0.547,0.632,0.707,1] ]    
    Mach = 20.
    shield = Newtonian(points,Mach)

    return shield
    
    
if __name__=="__main__":
    test=test_shield()
    test.analyse(5)
    #test.show()
    test.show_flight(0)
    print test.CA_T
    print test.CN_T
    print test.CM_T
    #test.CAa_plot()
    #test.CNa_plot()
    #test.CMa_plot()
    