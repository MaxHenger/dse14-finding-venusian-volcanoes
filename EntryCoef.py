# -*- coding: utf-8 -*-
"""
Created on Fri May 27 13:07:43 2016

This file contains the coefficients for the entry simulation
There are two variety, 2D and 3D.

@author: Julius
"""
import numpy as np
import Utility as util

class coef2D():
    def __init__(self):
        self.x=["V","y","R","p","q","r","a","b","m"]
        self.u=["de","da","dr","Mx","My","Mz"]
        
    def initial_size(self,m,S,b,c):
        self.m=m
        self.S=S
        self.b=b
        self.c=c
        
    def initial_MoI(self,Ixx,Iyy,Izz,Ixz):
        self.Ixx=Ixx
        self.Iyy=Iyy
        self.Izz=Izz
        self.Ixz=Ixz
        
        self.I=self.Ixx*self.Izz-self.Ixz**2
        self.Ip1=(self.Ixx-self.Iyy+self.Izz)/self.I
        self.Ip2=((self.Iyy-self.Izz)*self.Izz-self.Ixz**2)/self.I
        self.Ir1=((self.Ixx-self.Iyy)*self.Ixx-self.Ixz**2)/self.I
        self.Ir2=(-self.Ixx+self.Iyy-self.Izz)*self.Ixz/self.I
        
    
    def initial_misc(self,CL,L0,D0,g0):
        self.CL=CL
        self.D0=D0 # aerodynamic force
        self.L0=L0 # aerodynamic force
        self.g0=g0
        
    def initial_stab(self,CDa,CDM,CLM,CLa,Clb,CmM,Cma,Cnb,CSb,Clda,Cmda,Cnda,Cndr):
        self.CDa=CDa
        self.CDM=CDM
        self.CLM=CLM
        self.CLa=CLa
        self.Clb=Clb
        self.CmM=CmM
        self.Cma=Cma
        self.Cnb=Cnb
        self.CSb=CSb
        self.Clda=Clda
        self.Cmda=Cmda
        self.Cnda=Cnda
        self.Cndr=Cndr
        
    def initial(self,V0,y0,ydot0,R0,q0,p0,r0,a0,b0,m0):
        self.V0=V0
        self.y0=y0
        self.ydot0=ydot0
        self.R0=R0
        self.p0=p0
        self.q0=q0
        self.r0=r0
        self.a0=a0
        self.b0=b0
        self.m0=m0
        
        self.a=util.scale_a()
        self.M0=self.V0/self.a
        self.rho0=util.scale_height(R0)[2]
        self.q=0.5*self.rho0*self.V0**2 # dynamic pressure
        
    def matrixA(self):
        a= np.zeros((9,9))
        
        a[0][0]=self.aVV()
        a[0][1]=self.aVy()
        a[0][2]=self.aVR()
        a[0][3]=self.aVp()
        a[0][4]=self.aVq()
        a[0][5]=self.aVr()
        a[0][6]=self.aVa()
        a[0][7]=self.aVb()
        a[0][8]=self.aVm()
        
        a[1][0]=self.ayV()
        a[1][1]=self.ayy()
        a[1][2]=self.ayR()
        a[1][3]=self.ayp()
        a[1][4]=self.ayq()
        a[1][5]=self.ayr()
        a[1][6]=self.aya()
        a[1][7]=self.ayb()
        a[1][8]=self.aym()
        
        a[2][0]=self.aRV()
        a[2][1]=self.aRy()
        a[2][2]=self.aRR()
        a[2][3]=self.aRp()
        a[2][4]=self.aRq()
        a[2][5]=self.aRr()
        a[2][6]=self.aRa()
        a[2][7]=self.aRb()
        a[2][8]=self.aRm()
        
        a[3][0]=self.apV()
        a[3][1]=self.apy()
        a[3][2]=self.apR()
        a[3][3]=self.app()
        a[3][4]=self.apq()
        a[3][5]=self.apr()
        a[3][6]=self.apa()
        a[3][7]=self.apb()
        a[3][8]=self.apm()
        
        a[4][0]=self.aqV()
        a[4][1]=self.aqy()
        a[4][2]=self.aqR()
        a[4][3]=self.aqp()
        a[4][4]=self.aqq()
        a[4][5]=self.aqR()
        a[4][6]=self.aqa()
        a[4][7]=self.aqb()
        a[4][8]=self.aqm()
        
        a[5][0]=self.arV()
        a[5][1]=self.ary()
        a[5][2]=self.arR()
        a[5][3]=self.arp()
        a[5][4]=self.arq()
        a[5][5]=self.arr()
        a[5][6]=self.ara()
        a[5][7]=self.arb()
        a[5][8]=self.arm()
        
        a[6][0]=self.aaV()
        a[6][1]=self.aay()
        a[6][2]=self.aaR()
        a[6][3]=self.aap()
        a[6][4]=self.aaq()
        a[6][5]=self.aar()
        a[6][6]=self.aaa()
        a[6][7]=self.aab()
        a[6][8]=self.aam()
        
        a[7][0]=self.abV()
        a[7][1]=self.aby()
        a[7][2]=self.abr()
        a[7][3]=self.abp()
        a[7][4]=self.abq()
        a[7][5]=self.abr()
        a[7][6]=self.aba()
        a[7][7]=self.abb()
        a[7][8]=self.abm()
        
        a[8][0]=self.amV()
        a[8][1]=self.amy()
        a[8][2]=self.amR()
        a[8][3]=self.amp()
        a[8][4]=self.amq()
        a[8][5]=self.amr()
        a[8][6]=self.ama()
        a[8][7]=self.amb()
        a[8][8]=self.amm()
        
        self.MA=a
        return self.MA
    
    
#### Derivitives of stability
    # CD    
    def dCDdM(self):
        return self.CDM
    def dCDda(self):
        return self.CDa
    # CL 
    def dCLdM(self):
        return self.CLM
    def dCLda(self):
        return self.CLa
    def dCldb(self):
        return self.Clb
    # Cm    
    def dCmdM(self):
        return self.CmM
    def dCmda(self):
        return self.Cma
    # Cn    
    def dCndb(self):
        return self.Cnb
    # CS    
    def dCSdb(self):
        return self.CSb


#### functions cause I always type self.cos instead of np.cos
    def cos(self,angle):
        return np.cos(angle)
    def sin(self,angle):
        return np.sin(angle)
        
#### Derivatives with V
    def aVV(self):
        return -1/(self.m*self.V0)*(self.M0*self.dCDdM()*self.q*self.S+2*self.D0)
    def aVy(self):
        return -self.g0*np.cos(self.y0)
    def aVR(self):
        return 2* self.g0/self.R0*np.sin(self.y0)
    def aVp(self):
        return 0
    def aVq(self):
        return 0
    def aVr(self):
        return 0
    def aVa(self):
        return -1/self.m*self.dCDda()*self.q*self.S
    def aVb(self):
        return 0
    def aVm(self):
        return 0
        
#### Derivatives with y
    def ayV(self):
        return 1/self.V0*(-self.ydot0+2*self.V0/self.R0*np.cos(self.y0)) + np.cos(self.m0)/(self.m*self.V0**2)*(self.M0*self.dCLdM()*self.q*self.S+2*self.L0)
    def ayy(self):
        return -(self.V0/self.R0 - self.g0/self.V0)*np.sin(self.y0)
    def ayR(self):
        return (2*self.g0/self.V0-self.V0/self.R0)*np.cos(self.y0)/self.R0
    def ayp(self):
        return 0
    def ayq(self):
        return 0
    def ayr(self):
        return 0
    def aya(self):
        return np.cos(self.m0)/(self.m*self.V0)*self.dCLda()*self.q*self.S
    def ayb(self):
        return -np.sin(self.m0)/(self.m*self.V0)*self.dCSdb()*self.q*self.S
    def aym(self):
        return - self.L0/(self.m*self.V0)*np.sin(self.m0)
        
#### Derivatives with R
    def aRV(self):
        return np.sin(self.y0)
    def aRy(self):
        return self.V0*np.cos(self.y0)
    def aRR(self):
        return 0
    def aRp(self):
        return 0
    def aRq(self):
        return 0
    def aRr(self):
        return 0
    def aRa(self):
        return 0
    def aRb(self):
        return 0
    def aRm(self):
        return 0   
        
#### Derivatives with p
    def apV(self):
        return 0
    def apy(self):
        return 0
    def apR(self):
        return 0
    def app(self):
        return self.Ip2*self.q0
    def apq(self):
        return self.Ip1*self.p0 + self.Ip2*self.r0
    def apr(self):
        return self.Ip2*self.q0
    def apa(self):
        return 0
    def apb(self):
        return 1/self.Ixx *self.dCldb()*self.q*self.S*self.b
    def apm(self):
        return 0
        
#### Derivatives with q
    def aqV(self):
        return self.M0/(self.Iyy*self.V0)*self.dCmdM()*self.q*self.S*self.c
    def aqy(self):
        return 0
    def aqR(self):
        return 0
    def aqp(self):
        return -2*self.Ixz/self.Iyy*self.p0 + (self.Izz-self.Ixx)/self.Iyy*self.r0
    def aqq(self):
        return 0
    def aqr(self):
        return (self.Izz-self.Ixx)/self.Iyy*self.p0+2*self.Ixz/self.Iyy*self.r0
    def aqa(self):
        return 1/self.Iyy*self.dCmda()*self.q*self.S*self.c
    def aqb(self):
        return 0
    def aqm(self):
        return 0

#### Derivatives with r
    def arV(self):
        return 0
    def ary(self):
        return 0
    def arR(self):
        return 0
    def arp(self):
        return self.Ir1*self.q0
    def arq(self):
        return self.Ir1*self.p0+self.Ir2*self.r0
    def arr(self):
        return self.Ir2*self.q0
    def ara(self):
        return 0
    def arb(self):
        return 1/self.Izz*self.dCndb()*self.q*self.S*self.b
    def arm(self):
        return 0

#### Derivatives with a
    def aaV(self):
        return -self.g0/self.V0**2*np.cos(self.y0)*self.cos(self.m0)-1/(self.m*self.V0**2)*(self.M0*self.dCLdM()+self.CL)*self.q*self.S
    def aay(self):
        return -self.g0/self.V0*np.sin(self.y0)*np.cos(self.m0)
    def aaR(self):
        return 2*self.g0/(self.R0*self.V0)*np.cos(self.y0)*np.cos(self.m0)
    def aap(self):
        return 0
    def aaq(self):
        return 1
    def aar(self):
        return 0
    def aaa(self):
        return -1/(self.m*self.V0)*self.dCLda()*self.q*self.S
    def aab(self):
        return 0
    def aam(self):
        return -self.g0/self.V0*self.cos(self.y0)*self.sin(self.m0)
        
#### Derivatives with b
    def abV(self):
        return self.g0/self.V0**2*self.cos(self.y0)*self.sin(self.m0)
    def aby(self):
        return self.g0/self.V0*self.sin(self.y0)*self.sin(self.m0)
    def abR(self):
        return 2*self.g0/(self.R0*self.V0)*self.cos(self.y0)*self.sin(self.m0)
    def abp(self):
        return self.sin(self.a0)
    def abq(self):
        return 0
    def abr(self):
        return -self.cos(self.a0)
    def aba(self):
        return 0
    def abb(self):
        return -1/(self.m*self.V0)*self.dCSdb()*self.q*self.S
    def abm(self):
        return -self.g0/self.V0*self.cos(self.y0)*self.cos(self.m0)
        
#### Derivatives with m
    def amV(self):
        return np.tan(self.y0)*self.sin(self.m0)/(self.m*self.V0**2)*(self.M0*self.dCLdM()+self.CL)*self.q*self.S
    def amy(self):
        return self.L0/(self.m*self.V0)*self.sin(self.m0)
    def amR(self):
        return 0
    def amp(self):
        return -self.cos(self.a0)
    def amq(self):
        return 0
    def amr(self):
        return -self.sin(self.a0)
    def ama(self):
        return np.tan(self.y0)*self.sin(self.m0)/(self.m*self.V0)*self.dCLda()*self.q*self.S
    def amb(self):
        return np.tan(self.y0)*self.cos(self.m0)/(self.m*self.V0)*self.dCSdb()*self.q*self.S - self.L0/(self.m*self.V0)+self.g0/self.V0*self.cos(self.y0)*self.cos(self.m0)
    def amm(self):
        return np.tan(self.y0)*self.cos(self.m0)*self.L0/(self.m*self.V0)


    def matrixB(self):
        b=np.zeros((9,6))
        
        b[0][0]=self.bVe()
        b[0][1]=self.bVa()        
        b[0][2]=self.bVr()
        b[0][3]=self.bVx()
        b[0][4]=self.bVy()
        b[0][5]=self.bVz()
        
        b[1][0]=self.bye()
        b[1][1]=self.bya()        
        b[1][2]=self.byr()
        b[1][3]=self.byx()
        b[1][4]=self.byy()
        b[1][5]=self.byz()
        
        b[2][0]=self.bRe()
        b[2][1]=self.bRa()        
        b[2][2]=self.bRr()
        b[2][3]=self.bRx()
        b[2][4]=self.bRy()
        b[2][5]=self.bRz()
        
        b[3][0]=self.bpe()
        b[3][1]=self.bpa()        
        b[3][2]=self.bpr()
        b[3][3]=self.bpx()
        b[3][4]=self.bpy()
        b[3][5]=self.bpz()
        
        b[4][0]=self.bqe()
        b[4][1]=self.bqa()        
        b[4][2]=self.bqr()
        b[4][3]=self.bqx()
        b[4][4]=self.bqy()
        b[4][5]=self.bqz()
        
        b[5][0]=self.bre()
        b[5][1]=self.bra()        
        b[5][2]=self.brr()
        b[5][3]=self.brx()
        b[5][4]=self.bry()
        b[5][5]=self.brz()
        
        b[6][0]=self.bae()
        b[6][1]=self.baa()        
        b[6][2]=self.bar()
        b[6][3]=self.bax()
        b[6][4]=self.bay()
        b[6][5]=self.baz()
        
        b[7][0]=self.bbe()
        b[7][1]=self.bba()        
        b[7][2]=self.bbr()
        b[7][3]=self.bbx()
        b[7][4]=self.bby()
        b[7][5]=self.bbz()
        
        b[8][0]=self.bme()
        b[8][1]=self.bma()        
        b[8][2]=self.bmr()
        b[8][3]=self.bmx()
        b[8][4]=self.bmy()
        b[8][5]=self.bmz()
        
        self.MB=b
        return self.MB        
        
        
#### special derivitives
    def dCldda(self):
        return self.Clda
    def dCmdde(self):
        return self.Cmda
    def dCndda(self):
        return self.Cnda
    def dCnddr(self):
        return self.Cndr
        
        
#### Derivatives with V
    def bVe(self):
        return 0
    def bVa(self):
        return 0
    def bVr(self):
        return 0
    def bVx(self):
        return 0
    def bVy(self):
        return 0
    def bVz(self):
        return 0

#### Derivatives with y
    def bye(self):
        return 0
    def bya(self):
        return 0
    def byr(self):
        return 0
    def byx(self):
        return 0
    def byy(self):
        return 0
    def byz(self):
        return 0

#### Derivatives with R
    def bRe(self):
        return 0
    def bRa(self):
        return 0
    def bRr(self):
        return 0
    def bRx(self):
        return 0
    def bRy(self):
        return 0
    def bRz(self):
        return 0
        
#### Derivatives with p
    def bpe(self):
        return 0
    def bpa(self):
        return 1/self.Ixx*self.dCldda()*self.q*self.S*self.b
    def bpr(self):
        return 0
    def bpx(self):
        return self.Izz/self.I
    def bpy(self):
        return 0
    def bpz(self):
        return self.Ixz/self.I
        
#### Derivatives with q
    def bqe(self):
        return 1/self.Iyy*self.dCmdde()*self.q*self.S*self.c
    def bqa(self):
        return 0
    def bqr(self):
        return 0
    def bqx(self):
        return 0
    def bqy(self):
        return 1/self.Iyy
    def bqz(self):
        return 0
        
#### Derivatives with r
    def bre(self):
        return 0
    def bra(self):
        return 1/self.Izz*self.dCndda()*self.q*self.S*self.b
    def brr(self):
        return 1/self.Izz*self.dCnddr()*self.q*self.S*self.b
    def brx(self):
        return self.Ixz/self.I
    def bry(self):
        return 0
    def brz(self):
        return self.Ixx/self.I
        
#### Derivatives with a
    def bae(self):
        return 0
    def baa(self):
        return 0
    def bar(self):
        return 0
    def bax(self):
        return 0
    def bay(self):
        return 0
    def baz(self):
        return 0
        
#### Derivatives with b
    def bbe(self):
        return 0
    def bba(self):
        return 0
    def bbr(self):
        return 0
    def bbx(self):
        return 0
    def bby(self):
        return 0
    def bbz(self):
        return 0
        
#### Derivatives with m
    def bme(self):
        return 0
    def bma(self):
        return 0
    def bmr(self):
        return 0
    def bmx(self):
        return 0
    def bmy(self):
        return 0
    def bmz(self):
        return 0