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
    Vbal=mtot/(rhoatm-rhogas)*(buoyancyperc/100.)    
    mgas=rhogas*Vbal    
    return mtot,Vbal,mgas

def fullbuoyancy(mgas, mtot, Vbal, expancruise):    
    import VenusAtmosphere as atm
    #assume balloon goes down and balloonh was slightly expanded at higher alt. Adjest Volume for return to normal shape
    Vbalnew=Vbal*(100.-expancruise)/100.
    rhogasnew=mgas/Vbalnew

    #find what alt buoyancy is 100%
    #reiterate to find alt at which balloon lift is fully carrying mass
    found=False
    altbuoy=0
    while not found:
        Tatm, Patm, rhoatm, GravAcc=atm.VenusAtmosphere30latitude(altbuoy)        
        Lift=Vbalnew*(rhoatm-rhogasnew)        
        if (Lift-mtot)<1e-1:
            found=True
        else:
            altbuoy+=1e-2
    return (altbuoy)


def optimization_balloon_calc():
    mpayload = 90.
    hcruise = 50000
    m_molar = 2.016
    perc_buoy = 100
    n=50
    acc=15 # minimum 4
    import time
    start = time.time()
    for i in range(n):
        (ballooncalc(mpayload,hcruise,m_molar,perc_buoy))
    mid = time.time()
    for i in range(n):
        (math_ballooncalc(mpayload,hcruise,m_molar,perc_buoy,acc))
    end=time.time()
    print("Method A: ",mid-start)
    print("Method A: ",end-mid)
    print(ballooncalc(mpayload,hcruise,m_molar,perc_buoy))
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
    #random shit program
    
    if False:
        import VenusAtmosphere as atm
        hcruise = 50000 #km
        Cl =1.5
        mpayload = 90.
        V=40.
        hforce=hcruise/1000.
        mtotnew, Vbalnew,mgas=ballooncalc(mpayload, hforce, 2.016, 100)
        A=14
        print "here"
        found=False
        while not found:
            Tatm, Patm, rhoatm, GravAcc=atm.VenusAtmosphere30latitude(hforce)       
            rhogas = Patm/((8314.4598/2.016)*Tatm)
            Lift = GravAcc*Vbalnew*(rhoatm-rhogas)    
            altbuoy=fullbuoyancy(mgas,mtotnew,Vbalnew,10.)
            cord=(Vbalnew/2.4)**(1./3.)
            S=cord*A*cord
            LiftWing=1/2.*S*Cl*V**2*rhoatm
            if abs(Lift-LiftWing)<1:
                found=True
            else:
                hforce-=1e-2
        print hforce


    
    
