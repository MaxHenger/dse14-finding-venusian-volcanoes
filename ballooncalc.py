# -*- coding: utf-8 -*-
"""
Created on Tue May 03 11:29:28 2016

@author: Chaggai
"""


def ballooncalc(mpayload, hbal, molarmgas, buoyancyperc):
    import VenusAtmosphere as atm
    #get atmospheric contstants
    Tatm, Patm, rhoatm, GravAcc=atm.VenusAtmosphere30latitude(hbal)
    
    #calculate density volume and total mass of balloon 
    rhogas = Patm/((8314.4598/molarmgas)*Tatm)
    Vbal= mpayload/(rhoatm-rhogas)*(buoyancyperc/100.)    
    mtot=(rhogas*Vbal+mpayload)/0.2    
    
    #reiterate for snowball effect of balloon gas
    found=False
    while found==False:     
        Vbalnew=mtot/(rhoatm-rhogas)*(buoyancyperc/100.)
        mtotnew=(rhogas*Vbalnew+mpayload)/0.2        
        if abs(mtotnew-mtot)<1e-3:
            found=True
        else:
            mtot=mtotnew
    #calculate mass of gas in balloon 
    mgas=rhogas*Vbalnew    
    return(mtotnew, Vbalnew,mgas)
    
def math_ballooncalc(mpayload, hbal, molarmgas, buoyancyperc,accuracy=10):
    import VenusAtmosphere as atm
    Tatm, Patm, rhoatm, GravAcc=atm.VenusAtmosphere30latitude(hbal)
    #calculate density volume and total mass of balloon 
    rhogas = Patm/((8314.4598/molarmgas)*Tatm)
    mtot=sum([ mpayload*(1/0.2)**(i+1)*(buoyancyperc/100.)**(i)*(rhogas/(rhoatm-rhogas))**i for i in range(0,accuracy)]) 
    # iteration to determine fianl mass
    Vbal=mtot/(rhoatm-rhogas)*(buoyancyperc/100.)    
    mgas=rhogas*Vbal    
    return mtot,Vbal,mgas

def fullbuoyancy(mgas, mtot, Vbal, expancruise,stepSize=10,accuracy=10):    
    import VenusAtmosphere as atm
    import math
    #assume balloon goes down and balloon was slightly expanded at higher alt. Adjust Volume for return to normal shape
    Vbalnew=Vbal*(100.-expancruise)/100.
    rhogasnew=mgas/Vbalnew

    #find what alt buoyancy is 100%
    #Old Method: reiterate to find alt at which balloon lift is fully carrying mass
    #new Method: use scale height as estimation
    oldMethod=False
    if oldMethod:
        found=False
        altbuoy=0
        while not found:
            Tatm, Patm, rhoatm, GravAcc=atm.VenusAtmosphere30latitude(altbuoy)        
            Liftmass=Vbalnew*(rhoatm-rhogasnew)        
            if (Liftmass-mtot)<accuracy:
                found=True
            else:
                altbuoy+=stepSize
        return altbuoy
    else:
        rho0 = 65.
        H = 15.9*1000
        rhoatm=mtot/Vbalnew+rhogasnew 
        ScaleAlt= -math.log(rhoatm/rho0)*H
        
        altitude0 = ScaleAlt
        step0=100000   
        def findAlt(altitude,step):
            print(altitude,step)
            rhocurrent = atm.VenusAtmosphere30latitude(altitude)[2]
            if abs(rhocurrent-rhoatm)<accuracy:
                return altitude
            elif rhocurrent<rhoatm:
                return findAlt(altitude-step/2.,step/2.)
            else:
                return findAlt(altitude+step,step)
                
        return findAlt(altitude0,step0)
                
def optimization_balloon_calc(acc=15):
    """acc determines accuracy of method B, min 4, rec. 15 """
    mpayload = 90.
    hcruise = 50000
    m_molar = 2.016
    perc_buoy = 100
    n=50
    acc=acc # minimum 4
    import time
    start = time.time()
    for i in range(n):
        (ballooncalc(mpayload,hcruise,m_molar,perc_buoy))
    mid = time.time()
    for i in range(n):
        (math_ballooncalc(mpayload,hcruise,m_molar,perc_buoy,acc))
    end=time.time()
    print("Method A: ",mid-start)
    print(ballooncalc(mpayload,hcruise,m_molar,perc_buoy))
    print("Method A: ",end-mid)
    print(math_ballooncalc(mpayload,hcruise,m_molar,perc_buoy,acc))
    
    


def Solarpanelpower(Vbal, Thickcord, Aspect, alt):
    import Solar as sol
    #calc Apanel
    cord = (2.*Vbal/(Aspect*Thickcord))**(1./3.)
    span = Aspect*cord
    Ar=cord*span
    incmax=90.
    incmin=1.
    
    Psolarmax=sol.SolarPower(alt,Ar,incmax)
    Psolarmin=sol.SolarPower(alt,Ar,incmin)
    return(Psolarmax,Psolarmin)

if __name__=="__main__":
    
    mpayload = 90.
    hcruise = 50000
    m_molar = 2.016
    perc_buoy = 100
    acc=15 # minimum 4
    mtot,Vbal,mgas=math_ballooncalc(mpayload,hcruise,m_molar,perc_buoy,acc)
    print(fullbuoyancy(mgas,mtot,Vbal,10,accuracy=0.0001))
    
