# -*- coding: utf-8 -*-
"""
Created on Tue May 03 11:29:28 2016

@author: Chaggai
"""
import numpy as np
    
def balloonInital(mpayload=90, Hbouyancy=50000, molarmgas=2.016,structurePayloadFactor=0.2,accuracy=10):
    """ first order estimation of the balloon given bouyance altitude""" 
    import VenusAtmosphere as atm
    R = 8.314459848
    Tatm, Patm, rhoatm, GravAcc=atm.VenusAtmosphere30latitude(Hbouyancy)
    #calculate density volume and total mass of balloon 
    rhogas = Patm/((R/(molarmgas/1000.))*Tatm)
    mtot=sum([ mpayload*(1/structurePayloadFactor)**(i+1)*(rhogas/(rhoatm-rhogas))**i for i in range(0,accuracy)]) 
    # iteration to determine fianl mass
    Vbal=mtot/(rhoatm-rhogas)
    mgas=rhogas*Vbal    
    Presgas = rhogas*R/(molarmgas/1000.)*Tatm
    return mtot,molarmgas,Vbal,mgas,Presgas,Tatm

def balloonCruise(mtot,molarmgas,Vbal,mgas,Pgas,Tgas,expanFactor,contracFactor,cruiseBuoyancy,accuracy=0.0001):
    import VenusAtmosphere as atm
    R = 8.314459848
    Rsp = R/(molarmgas/1000.)
    mcruise=cruiseBuoyancy*mtot
    #rhocruise = P / (Rsp*T)
    altitude0=0
    step0=200000
    def findAlt(altitude,step):
        #assert altitude>0
        Tatm, Patm, rhoatm, GravAcc = atm.VenusAtmosphere30latitude(altitude)
        rhooutside = rhoatm

        #factor = sorted([1-contracFactor, Pgas/Patm * Tatm/Tgas, 1+expanFactor])[1]
        factor = max(min(1+expanFactor, Pgas/Patm * Tatm/Tgas), 1-contracFactor)
        
        Vcruise=factor*Vbal
        rhoinside=mgas/Vcruise
        
        deltaRho=mcruise/Vcruise
        currentDeltaRho=rhooutside-rhoinside
        #print(altitude,step,factor,deltaRho,currentDeltaRho)
        if abs(deltaRho-currentDeltaRho)<accuracy:
            return altitude,rhoinside
        elif deltaRho>currentDeltaRho:
            return findAlt(altitude-step/2.,step/2.)
        else:
            return findAlt(altitude+step/2.,step/2.)
            
    Hcruise,rhogas = findAlt(altitude0,step0)
    Tatm, Patm, rhoatm, GravAcc = atm.VenusAtmosphere30latitude(Hcruise)
    Pgas = rhogas*Rsp*Tatm
    return Hcruise,rhogas,Pgas,Patm
    

def Solarpanelpower(Sarea, alt,effectiveArea=0.7):
    import solar as sol
    #calc Apanel
    
    incmax=np.deg2rad(90.)
    incmin=np.deg2rad(1.)
    
    Psolarmax=sol.SolarPower(alt,Sarea,incmax)
    Psolarmin=sol.SolarPower(alt,Sarea,incmin)
    return(Psolarmax,Psolarmin)
    
def SolarMinInc(Sarea,Preq,alt,effectiveArea=0.7,accuracy=0.01):
    import solar as sol
    #calc Apanel
    angle0=0.
    step0=90.
    def PowerAvailiable(angle,step):
        inc=np.deg2rad(angle)
        Psolar=sol.SolarPower(alt,Sarea,inc)
        assert Preq-Psolar != nan
        if abs(Preq-Psolar)<accuracy:
            return inc
        elif Preq<Psolar:
            return PowerAvailiable(angle-step/2.,step/2.)
        else:
            return PowerAvailiable(angle+step/2.,step/2.)
        
    return np.rad2deg(PowerAvailiable(angle0,step0))
    
def PowerReqThrust(Thrust,velocity,area,density,powerFactor=0.15):
    powerIdeal =  0.5*Thrust*velocity*( ( Thrust/(area*velocity**2*density/2.) +1)**0.5 +1)
    #DV= Thrust*velocity
    #print(powerIdeal,DV)
    powerActual = powerIdeal*(1+powerFactor)
    
    return powerActual
    
    
def surfaceArea(Vbal,Thickcord,Aspect):
    cord = (2.*Vbal/(Aspect*Thickcord))**(1./3.)
    span = Aspect*cord
    Ar=cord*span #open area for folding and misc
    return Ar,cord,span
    
def SteadFlight(hcruise,cruiseBuoyancy,mtot,SurfaceArea,chord,Cl=0.5,Cd=0.08):
    #Lift = weight
    import VenusAtmosphere as atm
    Tatm, Patm, rhoatm, GravAcc = atm.VenusAtmosphere30latitude(hcruise)
    weight = mtot*abs(1-cruiseBuoyancy)*GravAcc
    velocity = (weight / (0.5*rhoatm*Cl*SurfaceArea) )**0.5
    drag = 0.5*rhoatm*velocity**2*SurfaceArea*Cd
    
    reynolds = velocity*rhoatm*chord/ DynViscocity(Tatm)

    return drag, velocity, rhoatm, reynolds
    
def DynViscocity(Temp,Viscinit=0.0000148,Tempinit=293.15,sutherland=240):
    a = 0.555*Tempinit+sutherland
    b = 0.555*Temp + sutherland
    mu=Viscinit*(a/b)*(Temp/Tempinit)**(3/2)
    return mu
    
if __name__=="__main__":
    
    mpayload = 100
    hbuoyancy = 50000
    structurePayloadFactor=0.2
    m_molar = 2.016
    
    expanFactor = 0.1
    contracFactor = 0.1
    cruiseBuoyancy=0.1
    
    Cl=0.5
    Cd=0.04
    ThickCord=0.14
    Aspect=12
    
    launcher_diameter=5.
    seperation=0.5
    n_engines=1.
    Pradius = (launcher_diameter - (1+n_engines)*seperation)/(2*n_engines)
    print(Pradius)
    Parea = n_engines*np.pi*Pradius**2
    Pfactor=0.15
    
    mtot,molarmgas,Vbal,mgas,pgas,tgas = balloonInital(mpayload,hbuoyancy,m_molar,structurePayloadFactor)
    Hcruise,rhogas_c,Pgas_c,Patm_c=balloonCruise(mtot,molarmgas,Vbal,mgas,pgas,tgas,expanFactor,contracFactor,cruiseBuoyancy)
    Sarea,chord,span = surfaceArea(Vbal,ThickCord,Aspect)
    Psolarmax,Psolarmin=Solarpanelpower(Sarea,Hcruise)
    drag,velocity,rho_c,reynolds= SteadFlight(Hcruise,cruiseBuoyancy,mtot,Sarea,chord,Cl,Cd)
    Ppower = PowerReqThrust(drag,velocity,Parea,rho_c,Pfactor)
    incmin = SolarMinInc(Sarea,Ppower,Hcruise)
    print(mtot,Hcruise,velocity,reynolds,Vbal,drag,Ppower,incmin)
