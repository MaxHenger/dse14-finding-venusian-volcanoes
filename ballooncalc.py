# -*- coding: utf-8 -*-
"""
Created on Tue May 03 11:29:28 2016

@author: Chaggai
"""
import numpy as np
    
#def math_ballooncalc_wrong(mpayload, hbal, molarmgas, buoyancyperc,accuracy=10):
#    """ first order estimation of the balloon """ 
#    import VenusAtmosphere as atm
#    R = 8.314459848
#    Tatm, Patm, rhoatm, GravAcc=atm.VenusAtmosphere30latitude(hbal)
#    #calculate density volume and total mass of balloon 
#    rhogas = Patm/((R/(molarmgas/1000.))*Tatm)
#    mtot=sum([ mpayload*(1/0.2)**(i+1)*(buoyancyperc/100.)**(i)*(rhogas/(rhoatm-rhogas))**i for i in range(0,accuracy)]) 
#    # iteration to determine fianl mass
#    Vbal=mtot/(rhoatm-rhogas)*(buoyancyperc/100.)    
#    mgas=rhogas*Vbal    
#    Presgas = rhogas*R/(molarmgas/1000.)*Tatm
#    return mtot,Vbal,mgas, Presgas
    
def balloonInital(mpayload=90, Hbouyancy=50000, molarmgas=2.016,structurPayloadFactor=0.2,accuracy=10):
    """ first order estimation of the balloon given bouyance altitude""" 
    import VenusAtmosphere as atm
    R = 8.314459848
    Tatm, Patm, rhoatm, GravAcc=atm.VenusAtmosphere30latitude(Hbouyancy)
    #calculate density volume and total mass of balloon 
    rhogas = Patm/((R/(molarmgas/1000.))*Tatm)
    mtot=sum([ mpayload*(1/structurPayloadFactor)**(i+1)*(rhogas/(rhoatm-rhogas))**i for i in range(0,accuracy)]) 
    # iteration to determine fianl mass
    Vbal=mtot/(rhoatm-rhogas)
    mgas=rhogas*Vbal    
    Presgas = rhogas*R/(molarmgas/1000.)*Tatm
    return mtot,molarmgas,Vbal,mgas,Presgas,Tatm

def balloonCruise(mtot,molarmgas,Vbal,mgas,Pgas,Tgas,expanFactor,contracFactor,cruiseBuoyancy,accuracy=0.00001):
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
    
    
def fullbuoyancy(mgas, mtot, Vbal, expancruise,accuracy=0.0001):    
    import VenusAtmosphere as atm
    #assume balloon goes down and balloon was slightly expanded at higher alt. Adjust Volume for return to normal shape
    Vbalnew=Vbal*(100.-expancruise)/100.
    rhogasnew=mgas/Vbalnew
    rhoatm=mtot/Vbalnew+rhogasnew 
    #find what alt buoyancy is 100%
    altitude0 = 0
    step0=200000   
    def findAlt(altitude,step):
        #print(altitude,step)
        rhocurrent = atm.VenusAtmosphere30latitude(altitude)[2]
        if abs(rhocurrent-rhoatm)<accuracy:
            return altitude
        elif rhocurrent<rhoatm:
            return findAlt(altitude-step/2.,step/2.)
        else:
            return findAlt(altitude+step/2.,step/2.)
            
    return findAlt(altitude0,step0)
                
    
def diveAltitude():
    pass

def Solarpanelpower(Vbal, Thickcord, Aspect, alt):
    import solar as sol
    #calc Apanel
    cord = (2.*Vbal/(Aspect*Thickcord))**(1./3.)
    span = Aspect*cord
    Ar=cord*span*0.7 #open area for folding and misc
    incmax=np.deg2rad(90.)
    incmin=np.deg2rad(1.)
    
    Psolarmax=sol.SolarPower(alt,Ar,incmax)
    Psolarmin=sol.SolarPower(alt,Ar,incmin)
    return(Psolarmax,Psolarmin,cord,span)
    


if __name__=="__main__":
    
    mpayload = 4500.
    hbuoyancy = 50000
    structurPayloadFactor=0.2
    m_molar = 2.016
    
    expanFactor = 0.1
    contracFactor = 0.1
    cruiseBuoyancy=0.1
    
    mtot,molarmgas,Vbal,mgas,pgas,tgas = balloonInital(mpayload,hbuoyancy,m_molar,structurPayloadFactor)
    print(mtot,Vbal,mgas)
    values=balloonCruise(mtot,molarmgas,Vbal,mgas,pgas,tgas,expanFactor,contracFactor,cruiseBuoyancy)
    #print(values)    
    #print(abs(values[3]-values[2]))
