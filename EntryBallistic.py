# -*- coding: utf-8 -*-
"""
Created on Mon Jun 06 11:22:16 2016

@author: Julius
"""

import numpy as np
import Utility as util
from matplotlib import pyplot as plt
import control
import Gravity 

import EntryNewtonian as mnf

class ballistic_sim:
    def __init__(self,shell,mass):
        self.shell=shell
        self.m=mass
        self.S=self.shell.SurfaceArea()
        self.c=2./3*self.shell.y[-1]
        self.diag=False
        self.rho_data=[]
        self.q_data=[]
        self.ydot_data=[]
        self.grav = Gravity.Gravity()
        
    def initial(self,V0,y0,R0,Iyy,q0,a0):
        self.para=False
        self.T_sim=0
        self.V0=V0
        self.y0=y0
        self.R0=R0
        self.q0=q0
        self.a0=a0
        self.Iyy=Iyy
        
        #self.update()
        
    def update(self):
        self.a=util.scale_a()
        self.M0=self.V0/self.a
        if self.diag:        
            print "V: ",self.V0
        if self.diag:        
            print "Mach: ", self.M0
        
        self.shell.update(self.M0)
        CT,CN,CM=shell.analyse(self.a0,mode="rad")
        self.CL=CN*np.cos(self.a0)-CT*np.sin(self.a0)
        self.CD=CT*np.cos(self.a0)+CN*np.sin(self.a0)
        if self.para and self.T_sim>=self.para_t:
            self.CD+=self.para_CD
            
        if self.diag:
            if self.para and self.T_sim>=self.para_t:
                print "Para Deployed: True"
            print "CL, CD: ",self.CL,self.CD
        
        self.q=0.5*util.scale_height(self.R0-Re)[2]*self.V0**2 # dynamic pressure
                
        if self.diag:        
            print "Dyn P: ",self.q
        
        if self.para and self.T_sim>=self.para_t:
            self.D0=self.q * self.S * (self.CD-self.para_CD) + self.q * self.para_S * self.para_CD
        else:
            self.D0=self.q * self.S * self.CD # aerodynamic force
        self.L0=self.q * self.S * self.CL # aerodynamic force
        if self.diag:        
            print "Lift, Drag : ",self.L0,self.D0
        
        self.g0=-self.grav(self.R0-Re,0,0)
        
        self.ydot0 = 0#-(self.V0/self.R0 - self.g0/self.V0)*np.cos(self.y0) + self.L0/(self.m*self.V0)
        
        #print self.ydot0        
        
        self.q_data.append(self.q)
        self.rho_data.append(util.scale_height(self.R0-Re)[2])
        self.ydot_data.append(self.ydot0)
        
        
    def MA(self,alpha,Mach):
        a = np.zeros((5,5))
        a[0][0]=self.aVV(alpha,Mach)
        a[0][1]=self.aVy()
        a[0][2]=self.aVR()
        a[0][3]=self.aVq()
        a[0][4]=self.aVa(alpha,Mach)
        
        a[1][0]=self.ayV(alpha,Mach)
        a[1][1]=self.ayy()
        a[1][2]=self.ayR()
        a[1][3]=self.ayq()
        a[1][4]=self.aya(alpha,Mach)
        
        a[2][0]=self.aRV()
        a[2][1]=self.aRy()
        a[2][2]=self.aRR()
        a[2][3]=self.aRq()
        a[2][4]=self.aRa()
        
        a[3][0]=self.aqV(alpha,Mach)
        a[3][1]=self.aqy()
        a[3][2]=self.aqR()
        a[3][3]=self.aqq()
        a[3][4]=self.aqa(alpha,Mach)
        
        a[4][0]=self.aaV(alpha,Mach)
        a[4][1]=self.aay()
        a[4][2]=self.aaR()
        a[4][3]=self.aaq()
        a[4][4]=self.aaa(alpha,Mach)
        
        return a
        
    def MB(self):
        return np.zeros((5,1))
        
    def dCDdM(self,alpha,Mach):
        return self.shell.CAM(alpha,Mach)
    def dCDda(self,alpha,Mach):
        return self.shell.CAa(alpha,Mach)
        
    def dCLdM(self,alpha,Mach):
        return self.shell.CNM(alpha,Mach)
    def dCLda(self,alpha,Mach):
        return self.shell.CNa(alpha,Mach)
        
    def dCmdM(self,alpha,Mach):
        return self.shell.CmM(alpha,Mach)
    def dCmda(self,alpha,Mach):
        return self.shell.Cma(alpha,Mach)
    
    
    
#### Derivatives with V
    def aVV(self,alpha,Mach):
        return -1./(self.m*self.V0)*(self.M0*self.dCDdM(alpha,Mach)*self.q*self.S+2*self.D0)
    def aVy(self):
        return -self.g0*np.cos(self.y0)
    def aVR(self):
        return 2* self.g0/self.R0*np.sin(self.y0)  
    def aVq(self):
        return 0
    def aVa(self,alpha,Mach):
        return -1./self.m*self.dCDda(alpha,Mach)*self.q*self.S
        
#### Derivatives with y
    def ayV(self,alpha,Mach):
        return 1./self.V0*(-self.ydot0+2*self.V0/self.R0*np.cos(self.y0)) + 1./(self.m*self.V0**2)*(self.M0*self.dCLdM(alpha,Mach)*self.q*self.S+2*self.L0)
    def ayy(self):
        return 0#0.1*-(self.V0/self.R0 - self.g0/self.V0)*np.sin(self.y0)
    def ayR(self):
        return (2*self.g0/self.V0-self.V0/self.R0)*np.cos(self.y0)/self.R0
    def ayq(self):
        return 0
    def aya(self,alpha,Mach):
        return 1./(self.m*self.V0)*self.dCLda(alpha,Mach)*self.q*self.S        
        
#### Derivatives with R
    def aRV(self):
        return np.sin(self.y0)
    def aRy(self):
        return self.V0*np.cos(self.y0)
    def aRR(self):
        return 0
    def aRq(self):
        return 0
    def aRa(self):
        return 0
        
#### Derivatives with q
    def aqV(self,alpha,Mach):
        if self.diag:
            print "alpha: ",alpha
            print "Mach: ",Mach
            print "dCmdM: ",self.dCmdM(alpha,Mach)
            print "TEST: ",self.M0/(self.Iyy*self.V0)*self.dCmdM(alpha,Mach)*self.q*self.S*self.c
        return self.M0/(self.Iyy*self.V0)*self.dCmdM(alpha,Mach)*self.q*self.S*self.c
    def aqy(self):
        return 0
    def aqR(self):
        return 0
    def aqq(self):
        return 0
    def aqa(self,alpha,Mach):
        if self.diag:       
            print "alpha: ",alpha
            print "Mach: ",Mach
            print "dCmda: ",self.dCmda(alpha,Mach)
            print "TEST: ",1./self.Iyy*self.dCmda(alpha,Mach)*self.q*self.S*self.c
        return 1./self.Iyy*self.dCmda(alpha,Mach)*self.q*self.S*self.c
        
#### Derivatives with a
    def aaV(self,alpha,Mach):
        return -self.g0/self.V0**2*np.cos(self.y0)-1./(self.m*self.V0**2)*(self.M0*self.dCLdM(alpha,Mach)+self.CL)*self.q*self.S
    def aay(self):
        return -self.g0/self.V0*np.sin(self.y0)
    def aaR(self):
        return -2*self.g0/(self.R0*self.V0)*np.cos(self.y0)
    def aaq(self):
        return 1.
    def aaa(self,alpha,Mach):
        return -1./(self.m*self.V0)*self.dCLda(alpha,Mach)*self.q*self.S
        
    def simulate(self,length=50,len_block=0.1,dt=0.1,restart=False):
        self.dt=dt
        
        if restart:
            print "Restarting Simulation\n"
            self.youts=[]
            self.T = np.arange(0,length,dt)
        elif "youts" not in dir(self):
            print "Previous results not found\nRestarting simulatios\n"
            self.youts=[]
            self.T = np.arange(0,length,dt)
        else:
            print "Continuing Simulation"
            self.T = np.hstack((self.T,np.arange(self.T[-1]+dt,self.T[-1]+dt+length,dt)))
            
        count_blocks=int(length/len_block)
        print "Amount of blocks: ",count_blocks
        
        for i in range(0,count_blocks):
            T0=i*len_block
            if (i+1)%50==0:
                print "\nDoing Block: ",i+1,"/",count_blocks
            if len(self.youts)!=0:
                init=self.youts[-1][-1].T # [V0,y0,R0,q0,a0]
            else:
                init=[self.V0,self.y0,self.R0,self.q0,self.a0]
            try:
                self.youts.append(self.simulate_block(init,T0,len_block,dt))
            except:
                print("Failed at block: ",i+1)
                raise
        
    def simulate_block(self,init,T0,len_block,dt):
        T = np.arange(T0-dt,T0+len_block,dt)
        u = np.zeros((len(T)))
        if self.diag:        
            print "Time segment: "+str(T[0])+" to "+str(T[-1])
        self.V0,self.y0,self.R0,self.q0,self.a0 = init
        self.update()
        MA = self.MA(self.a0,self.M0)
        MB = self.MB()
        MC = np.eye(len(MA))
        MD = np.zeros(MB.shape)        
        ssS=control.ss(MA,MB,MC,MD)
        if self.diag:        
            print MA
            print ssS
        yout,T,xout = control.lsim(ssS,u,T,init)     
        fix1=yout[1:]
        self.T_sim=T[-1]
        return fix1[:int(len_block/dt)]
        
    def calc_accel(self,V,dt=0.1):
        self.accel=[]
        for i in range(1,len(V)-1):
            self.accel.append( (V[i+1]-V[i-1])/(2*dt) )
        return self.accel
        
    def plot(self,yout=None,T=None,Rtrim=0,Ltrim=0):
        if yout==None:
            y=self.youts[0]
            for i in range(1,len(self.youts)):
                y=np.vstack((y,self.youts[i]))
            T = self.T
        else: 
            y = yout.T
        y = y.T
        assert Ltrim%1==0
        assert Rtrim%1==0
        if Ltrim==0:
            y = y[:,Rtrim:]
        else:
            y = y[:,Rtrim:-Ltrim]
        print "Length Time vs Data", len(T)," vs ",len(y[0])
        if len(T)>len(y[0]):
            print "Adjusting Time span"
            T = np.arange(0,len(y[0])*self.dt,self.dt)
        
        figd=plt.figure("Entry Simulation Results", figsize=(20,11))
        ncols=2
        nrows=3
        n=0
        ax=[]
        for i in range(nrows):
                for j in range(ncols):
                    n+=1
                    ax.append( figd.add_subplot(nrows*100+ncols*10+n) )
            
        n=0
        plt.rcParams.update({"font.size":14}) 
        ax[n].plot(T, y[n], label = r"Velocity", color="blue", linestyle="-")
        #ax[n].plot(T, roll_sim,label = r"Numerical improved", color="red", linestyle= "--" )
        #ax[n].plot(T, roll_sim_old,label = r"Numerical original", color="green",linestyle=":", linewidth=2.5)
        ax[n].set_title(r"Velocity")
        ax[n].set_ylabel(r"Velocity in $m/s$", fontsize=16)
        
        n=1
        ax[n].plot(T, y[n]*180./np.pi, label = r"flight path", color="blue", linestyle="-")
        #ax[n].plot(T, yaw_sim,label = r"Numerical improved", color="red", linestyle= "--" )
        #ax[n].plot(T, yaw_sim_old,label = r"Numerical original", color="green",linestyle=":", linewidth=2.5)
        ax[n].set_title(r"Flight path angle")
        ax[n].set_ylabel(r"$\gamma$ in $deg$", fontsize=16)
        n=2
        ax[n].plot(T, y[2]-Re, label = r"Altitude", color="blue", linestyle="-")
        #ax[n].plot(T, rollrate_sim,label = r"Numerical improved", color="red", linestyle= "--" )
        #ax[n].plot(T, rollrate_sim_old,label = r"Numerical original", color="green",linestyle=":", linewidth=2.5)
        ax[n].set_title(r"Altitude")
        ax[n].set_ylabel(r"$h$ in $m$", fontsize=16)
        n=3
        ax[n].plot(T, y[-2]*180./np.pi, label = r"pitch rate", color="blue", linestyle="-")
        #ax[n].plot(T, rollrate_sim,label = r"Numerical improved", color="red", linestyle= "--" )
        #ax[n].plot(T, rollrate_sim_old,label = r"Numerical original", color="green",linestyle=":", linewidth=2.5)
        ax[n].set_title(r"pitch rate")
        ax[n].set_ylabel(r"$q$ in $deg/s$", fontsize=16)
        n=4
        ax[n].plot(T[:-2], self.calc_accel(y[0],self.dt), label = r"acceleration", color="blue", linestyle="-")
        #ax[n].plot(T, rollrate_sim,label = r"Numerical improved", color="red", linestyle= "--" )
        #ax[n].plot(T, rollrate_sim_old,label = r"Numerical original", color="green",linestyle=":", linewidth=2.5)
        ax[n].set_title(r"Acceleration")
        ax[n].set_ylabel(r"$a$ in $m/s^2$", fontsize=16)
        n=5
        ax[n].plot(T, y[-1]*180./np.pi, label = r"pitch", color="blue", linestyle="-")
        #ax[n].plot(T, rollrate_sim,label = r"Numerical improved", color="red", linestyle= "--" )
        #ax[n].plot(T, rollrate_sim_old,label = r"Numerical original", color="green",linestyle=":", linewidth=2.5)
        ax[n].set_title(r"pitch")
        ax[n].set_ylabel(r"$\alpha$ in $deg$", fontsize=16)
        
        print "Max g: ",np.min(self.calc_accel(y[0],self.dt))/9.81
    
    def parachute(self,Time,S,CD):
        self.para=True
        self.para_t=Time
        self.para_S=S
        self.para_CD=CD
        
    def _import(self,FilePath=".\sim.save"):
        import pickle
        self.youts=pickle.load( open( FilePath, "rb" ) )
        
    def _export(self,FilePath=".\sim.save"):
        import pickle
        pickle.dump(self.youts,open(FilePath,"wb"))
        
        
if __name__=="__main__":
    shell = mnf.test_shield()
    
    mass = 1500.
    
    y0      = -10 *np.pi/180     # flight path angle [rad]
    V0      = 7260.            # Inital Velocity [m/s]
    Re      = 6052.*1000        # Radius of venus [m]
    R0      = Re + 260*1000     # Inital Orbital Radius [m]
    
    q0      =   0               # pitch rate [rad/s]
    a0      =   0               # pitch angle [rad]
    
    Iyy     = 6000.           # kg / m^2
        
    
    sim = ballistic_sim(shell,mass)
    sim.initial(V0,y0,R0,Iyy,q0,a0) # order: V0,y0,ydot0,R0,Iyy,q0,a0
    sim.diag=False
    sim.parachute(57,50,0.8)
    sim.simulate(length=300,len_block=0.1,dt=0.1,restart=True)
    #sim.simulate(length=5,len_block=0.1,dt=0.01,restart=False)
    sim.plot()