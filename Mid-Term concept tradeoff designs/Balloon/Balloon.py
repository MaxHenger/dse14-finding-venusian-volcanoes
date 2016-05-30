# -*- coding: utf-8 -*-
"""
Created on Mon May 09 16:00:58 2016

@author: Dong-ho
"""

import numpy as np
import os

os.chdir('../../')
    
def balloonInital(mpayload, Hbouyancy, molarmgas,structurPayloadFactor,accuracy):
    """ first order estimation of the balloon given bouyance altitude""" 
    import VenusAtmosphere as atm
    R = 8.314459848
    Tatm, Patm, rhoatm, GravAcc=atm.VenusAtmosphere30latitude(Hbouyancy)
    #calculate density volume and total mass of balloon 
    rhogas = Patm/((R/(molarmgas/1000.))*Tatm)
    mtot=169.
    #sum([ mpayload*(1/structurPayloadFactor)**(i+1)*(rhogas/(rhoatm-rhogas))**i for i in range(0,accuracy)]) 
    # iteration to determine final mass
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
    
    
#def fullbuoyancy(mgas, mtot, Vbal, expancruise,accuracy=0.0001):    
#    import VenusAtmosphere as atm
#    #assume balloon goes down and balloon was slightly expanded at higher alt. Adjust Volume for return to normal shape
#    Vbalnew=Vbal*(100.-expancruise)/100.
#    rhogasnew=mgas/Vbalnew
#    rhoatm=mtot/Vbalnew+rhogasnew 
#    #find what alt buoyancy is 100%
#    altitude0 = 0
#    step0=200000   
#    def findAlt(altitude,step):
#        #print(altitude,step)
#        rhocurrent = atm.VenusAtmosphere30latitude(altitude)[2]
#        if abs(rhocurrent-rhoatm)<accuracy:
#            return altitude
#        elif rhocurrent<rhoatm:
#            return findAlt(altitude-step/2.,step/2.)
#        else:
#            return findAlt(altitude+step/2.,step/2.)
#            
#    return findAlt(altitude0,step0)
                

def Solarpanelpower(Vbal, alt):
    import solar as sol
    #calc Apanel
#    cord = (2.*Vbal/(Aspect*Thickcord))**(1./3.)
#    span = Aspect*cord
#    Ar=cord*span*0.7 #open area for folding and misc

    r =  (3./(4.*np.pi)*Vbal)**(1./3.)  
    Shalf = 2.*np.pi*r**2
    
    incmax=np.deg2rad(90.)
    incmin=np.deg2rad(1.)
    
    Psolarmax=sol.SolarPower(alt,Shalf,incmax)
    Psolarmin=sol.SolarPower(alt,Shalf,incmin)
    return(Psolarmax,Psolarmin,r,Shalf)
    


if __name__=="__main__":
    
    mpayload = 15.
    structurePayloadFactor = 15./(159.-15.-16.)
    accuracy = 50
#    hbuoyancy = 50000
  
    m_molar_H2 = 2.01588
    m_molar_He = 4.002602
    m_molar_NH3 = 17.031
    m_molar_H20 = 18.01528
    m_molar_N2 = 28.0134
    
    m_molar_tab = [m_molar_H2, m_molar_He, m_molar_NH3, m_molar_H20, m_molar_N2]
    alt_tab = [i for i in range(0,105000,5000)] 
    
    specs = []
    for i in range(len(alt_tab)):
        specs.append([alt_tab[i]])
        
    for i in range(len(alt_tab)):
        for j in range(len(m_molar_tab)):
            mtot,molarmgas,Vbal,mgas,pgas,tgas = balloonInital(mpayload,alt_tab[i],m_molar_tab[j],structurePayloadFactor,accuracy)
            Psolarmax,Psolarmin,r,Shalf = Solarpanelpower(Vbal,alt_tab[i])
            
            specs[i].append([r,Vbal,mtot,Shalf,Psolarmax,Psolarmin])
            
    for i in range(len(alt_tab)):
        print 'ALTITUDE = ',specs[i][0]
        print '------------------------------------------------------------ '
        
        for j in range(len(m_molar_tab)):
            print '---GAS ',j+1, '---'
            print 'Radius = ',specs[i][j+1][0], 'm'
            print 'Volume = ',specs[i][j+1][1], 'm^3'
            print 'Mass = ',specs[i][j+1][2], 'kg'
            print 'Solar panel area = ',specs[i][j+1][3], 'm^2'
            print 'P solar max = ',specs[i][j+1][4], 'W'
            print 'P solar min = ',specs[i][j+1][5], 'W'
            print ' '

#    print 'ALTITUDE = ',specs[11][0]
#    print '------------------------------------------------------------ '
#    j=1    
#    print '---GAS ',j+1, '---'
#    print 'Radius = ',specs[i][j+1][0], 'm'
#    print 'Volume = ',specs[i][j+1][1], 'm^3'
#    print 'Mass = ',specs[i][j+1][2], 'kg'
#    print 'Solar panel area = ',specs[i][j+1][3], 'm^2'
#    print 'P solar max = ',specs[i][j+1][4], 'W'
#    print 'P solar min = ',specs[i][j+1][5], 'W'
#    print ' '    
    
#    structurePayloadFactor=0.2
#    m_molar = 2.016
#    
#    expanFactor = 0.1
#    contracFactor = 0.1
#    cruiseBuoyancy=0.1
#    
#    Cl=0.5
#    Cd=0.08
#    ThickCord=0.14
#    Aspect=12
#    
#    Pradius = 0.5
#    n_engines=4
#    Parea = n_engines*np.pi*Pradius**2
#    Pfactor=0.15
    
#    mtot,molarmgas,Vbal,mgas,pgas,tgas = balloonInital(mpayload,hbuoyancy,m_molar,structurePayloadFactor)
#    Hcruise,rhogas_c,Pgas_c,Patm_c=balloonCruise(mtot,molarmgas,Vbal,mgas,pgas,tgas,expanFactor,contracFactor,cruiseBuoyancy)
#    Psolarmax,Psolarmin=Solarpanelpower(Vbal,Hcruise)
#    Sarea,cord,span = surfaceArea(Vbal,ThickCord,Aspect)
#    drag,velocity,rho_c= SteadFlight(Hcruise,cruiseBuoyancy,mtot,Sarea,Cl,Cd)
#    Ppower = PowerReqThrust(drag,velocity,Parea,rho_c,Pfactor)
#    print(mtot,Hcruise,Psolarmin,Psolarmax)
    
    

    
    