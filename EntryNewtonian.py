# -*- coding: utf-8 -*-
"""
Created on Thu Jun 02 16:23:00 2016

@author: Julius
"""

import numpy as np
import matplotlib.pyplot as plt

class Newtonian():
    def __init__(self,points):
        self.geometry(points)
        self.gam = 1.2941
        
    def update(self,Mach):
        self.Mach=Mach
        limit=1.1
        if Mach<=limit:
            self.Mach=limit
        try:
            assert self.Mach>=1
        except AssertionError:
            print "Mach is: ",Mach
            raise
        self.CP_max = self.CP(self.Mach,self.gam)
        
    def SurfaceArea(self):
        return self.y[-1]**2*np.pi
        
    def CP(self,Mach,gamma):
        return ( (self.pfrac(Mach,gamma)) -1 )/( (gamma/2.)*Mach**2)
        #return (((1+gamma)**2/(4*gamma))**(gamma/(gamma-1)) * 4/(gamma+1))
        #return 2./(gamma*Mach**2)*( ( ((gamma+1)**2*Mach**2)/(4*gamma*Mach**2-2*(gamma-1))**(gamma/(gamma-1)) ) * (1-gamma+2*gamma*Mach**2)/(gamma+1) - 1  )
    
    def pfrac(self,Mach,gamma):
        #return (2*gamma*Mach**2 - (gamma-1))/(gamma+1)
        # https://books.google.nl/books?id=wGsTBwAAQBAJ&pg=PA70&lpg=PA70&dq=lees+modified+newtonian+flow&source=bl&ots=PSB1WZEJeD&sig=bi4iZEoG57AFOkuhTdSWQGQw-J0&hl=en&sa=X&redir_esc=y#v=onepage&q=lees%20modified%20newtonian%20flow&f=false
        return ((1+gamma)**2 *Mach**2 /(4*gamma*Mach**2 - 2*(gamma-1)))**(gamma/(gamma-1)) * ((1-gamma+2*gamma*Mach**2)/(gamma+1))
        
    def gamma(self,rho0,rho1,Mach):
        #NOT used
        n=rho1/rho0
        self.gam = (n+1)/(n-1) - 2*n/( (n-1)**Mach**2 )
        
    def CA(self,theta,alpha):
        #return 2*np.sin(theta)**2 + np.sin(alpha)**2*(1-3*np.sin(theta)**2)
        CA = self.CP_max*np.sin(theta)**2 + np.sin(alpha)**2*(1-3*np.sin(theta)**2)
        if alpha>theta:
            #print "CA: alpha > theta"
            beta = np.arcsin(np.tan(theta)/np.tan(alpha))
            cosb = np.cos(beta)
            T2 = (beta+np.pi/2.0)/np.pi
            T4 = cosb * np.sin(2*alpha) * np.sin(2*theta)
            CA = CA*T2 + 3*T4/(4*np.pi)
            
        return CA
        
    def CN(self,theta,alpha):
        #CN =  np.cos(theta)**2*np.sin(2*alpha)
        CN = self.CP_max/2.*np.cos(theta)**2*np.sin(2*alpha)
        if alpha>theta:
            #print "CN: alpha > theta"
            beta = np.arcsin(np.tan(theta)/np.tan(alpha))
            cosb = np.cos(beta)
            T2 = (beta+np.pi/2.0)/np.pi
            T3 = np.tan(theta)/np.tan(alpha)
            T5 = cosb/(3.*np.pi)
            CN = CN*(T2+T5*(T3+2./T3))
        return CN
        
    def CM(self,theta,alpha):
        return -2./3*(self.CN(theta,alpha)/(np.tan(theta)*np.cos(theta)**2 ))
        #return -self.CP_max/3.*(self.CN(theta,alpha)/(np.tan(theta)*np.cos(theta)**2 ))
        
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
        
    def CAM(self,alpha,Mach,dM=0.1):
        self.update(Mach-dM)
        CA1,CN1,CM1 = self.analyse(alpha,mode="rad")
        CD1 = self.CD(CA1,CN1,alpha)
        self.update(Mach+dM)
        CA2,CN2,CM2 = self.analyse(alpha,mode="rad")
        CD2 = self.CD(CA2,CN2,alpha)
        return (CD2-CD1)/(2*dM)
        
    def CAa(self,alpha,Mach,da=0.01):
        self.update(Mach)
        CA1,CN1,CM1 = self.analyse(alpha,mode="rad")
        CA2,CN2,CM2 = self.analyse(alpha,mode="rad")
        CD1 = self.CD(CA1,CN1,alpha)
        CD2 = self.CD(CA2,CN2,alpha)
        return (CD2-CD1)/(2*da)
        
    def CNM(self,alpha,Mach,dM=0.1):
        self.update(Mach-dM)
        CA1,CN1,CM1 = self.analyse(alpha,mode="rad")
        CL1 = self.CL(CA1,CN1,alpha)
        self.update(Mach+dM)
        CA2,CN2,CM2 = self.analyse(alpha,mode="rad")
        CL2 = self.CL(CA2,CN2,alpha)
        return (CL2-CL1)/(2*dM)
    
    def CNa(self,alpha,Mach,da=0.01):
        self.update(Mach)
        CA1,CN1,CM1 = self.analyse(alpha,mode="rad")
        CA2,CN2,CM2 = self.analyse(alpha,mode="rad")
        CL1 = self.CL(CA1,CN1,alpha)
        CL2 = self.CL(CA2,CN2,alpha)
        return (CL2-CL1)/(2*da)
        
    def CmM(self,alpha,Mach,dM=0.01):
        self.update(Mach-dM)
        Cm1 = self.analyse(alpha,mode="rad")[2]
        self.update(Mach+dM)
        Cm2 = self.analyse(alpha,mode="rad")[2]
        return (Cm2-Cm1)/(2*dM)
    
    def Cma(self,alpha,Mach,da=0.01):
        self.update(Mach)
        Cm1 = self.analyse(alpha-da,mode="rad")[2]
        Cm2 = self.analyse(alpha+da,mode="rad")[2]
        return (Cm2-Cm1)/(2*da)
                
    def geometry(self,points):
        self.points=np.array(points)
        self.x=self.points[0]
        self.y=self.points[1]
        self.dx=self.x[1:]-self.x[:-1]
        self.dy=self.y[1:]-self.y[:-1]
        self.theta=np.arctan( self.dy/self.dx )
        self.s = self.y[1:]/(self.dy/self.dx)
        
    
        
        
##### PLOTTING FUNCTIONS #####        
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
        
    def CAa_plot(self,Mach=15):
        self.update(Mach)
        a = np.arange(-5,10,0.1)
        CAa=[]
        for i in a:
            CAa.append(self.analyse(i)[0])
        plt.plot(a,CAa)
        plt.ylabel(r"$C_A [-]$",fontsize=14)
        plt.xlabel(r"$\alpha$ [deg]",fontsize=14)
        plt.grid(True)
        plt.show()
    
    def CAa_plot_report(self):
        M = [2,5,10,25]
        linestyles=["--",":","-.","-"]
        linewidths =[1,2,1.3,1]
        plt.figure(figsize=(5,4))
        for i,M_it in enumerate(M):
            self.update(M_it)
            a = np.arange(-5,40,0.1)
            CAa=[]
            for a_it in a:
                CAa.append(self.analyse(a_it)[0])
            plt.plot(a,CAa,linewidth=linewidths[i],linestyle=linestyles[i],label="Mach: "+str(M_it))
            
        plt.ylabel(r"$C_A \; [-]$", fontsize=14)
        plt.xlabel(r"$\alpha \; [deg]$",fontsize=14)
        plt.grid(True)
        plt.legend(loc=3)
        plt.tight_layout()
        plt.show()
        
    def CNa_plot(self,Mach=15):
        self.update(Mach)
        a = np.arange(-5,10,0.1)
        CNa=[]
        for i in a:
            CNa.append(self.analyse(i)[1])
        plt.plot(a,CNa)
        plt.show()
        
    def CNa_plot_report(self):
        M = [2,5,10,25]
        linestyles=["--",":","-.","-"]
        linewidths =[1,2,1.3,1]
        plt.figure(figsize=(5,4))
        for i,M_it in enumerate(M):
            self.update(M_it)
            a = np.arange(-5,40,0.1)
            CNa=[]
            for a_it in a:
                CNa.append(self.analyse(a_it)[1])
            plt.plot(a,CNa,linewidth=linewidths[i],linestyle=linestyles[i],label="Mach: "+str(M_it))
            
        plt.ylabel(r"$C_N \; [-]$", fontsize=14)
        plt.xlabel(r"$\alpha \; [deg]$",fontsize=14)
        plt.grid(True)
        plt.legend(loc=2)
        plt.tight_layout()
        plt.show()
        
    def CMa_plot(self,Mach=15):
        self.update(Mach)
        a = np.arange(-10,50,0.1)
        CMa=[]
        for i in a:
            CMa.append(self.analyse(i)[2])
        plt.plot(a,CMa)
        plt.show()
        
    def CMa_plot_report(self):
        M = [2,5,10,25]
        linestyles=["--",":","-.","-"]
        linewidths =[1,2,1.3,1]
        plt.figure(figsize=(5,4))
        for i,M_it in enumerate(M):
            self.update(M_it)
            a = np.arange(-5,40,0.1)
            CMa=[]
            for a_it in a:
                CMa.append(self.analyse(a_it)[2])
            plt.plot(a,CMa,linewidth=linewidths[i],linestyle=linestyles[i],label="Mach: "+str(M_it))
            
        plt.ylabel(r"$C_M \; [-]$", fontsize=14)
        plt.xlabel(r"$\alpha \; [deg]$",fontsize=14)
        plt.grid(True)
        plt.legend(loc=1)
        plt.tight_layout()
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
        M = np.arange(0.5,20,0.5)
        CM=[]
        for i in M:
            self.CP_max = self.CP(i,self.gam)
            CM.append(self.analyse(alpha)[2])
        self.CP_max = self.CP(self.Mach,self.gam)
        plt.plot(M,CM)
        plt.show()
        
    def CP_plot(self):
        M = np.arange(0,20,0.5)
        CP=[]
        for i in M:
            CP.append(self.CP(i,self.gam))
        self.CP_max = self.CP(self.Mach,self.gam)
        plt.plot(M,CP)
        plt.show()
        
    def CL(self,CA,CN,alpha):
        return CN*np.cos(alpha)-CA*np.sin(alpha)        
        
    def CD(self,CA,CN,alpha):
        return CN*np.cos(alpha)+CN*np.sin(alpha)
        
        
def test_shield():
    depth = 1.5
    width = 4.6
    depth_nose = 0.05
    width_nose = 1
    
    import math
    RN = 1.15
    X = []
    Y = []
    for x in np.arange(0, 1.82,0.01):
        if x <= RN - RN * math.cos(0.25*math.pi):
            y = RN * math.sin(np.arccos((RN-x)/RN))
            X.append(x)
            Y.append(y)
        elif x > RN - RN * math.cos(0.25*math.pi):
            y = x + RN*math.sin(0.25*math.pi) - (RN - RN * math.cos(0.25*math.pi))
            X.append(x)
            Y.append(y)
    
    ##plot shape
#    plt.plot(X, Y)
#    plt.axis([0, X[-1], 0, Y[-1]])
#    plt.show()
#    def y(x):
#        if x <= depth_nose:
#            return width_nose/2.*(x/depth_nose)**0.3
#        else:
#            return width_nose/2.+width/2.*((x-depth_nose)/depth)**0.8
#    dt=0.01
#    x = np.arange(0,depth+dt,dt)
#    yout=np.zeros(len(x))
#    for i,x_i in enumerate(x):
#        yout[i]=y(x_i)
#    points=np.array([x,yout])
    #print x
    #print yout
    Mach=10
    #points=[[0,0.1,0.2,0.3,0.4,0.5,2],[0,0.3162,0.447,0.547,0.632,0.707,1] ]    
    shield = Newtonian([X,Y])
    shield.update(Mach)
    return shield
    
    
if __name__=="__main__":
    test=test_shield()
    #test.analyse(2)
    #test.show()
    #test.show_flight(0)
    #print test.CA_T
    #print test.CN_T
    #print test.CM_T
    #test.CAa_plot()
    #test.CNa_plot()
    #test.CMa_plot(10)
    test.CAa_plot_report()
    