# -*- coding: utf-8 -*-
"""
Created on Tue May 03 11:29:28 2016

@author: Chaggai
"""


def ballooncalc(mpayload, hbal, molarmgas, buoyancyperc):
    import VenusAtmosphere as atm
    #get atmospheric contstants
    Tatm, Patm, rhoatm, GravAcc=atm.VenusAtmosphere30latitude(hbal)
    Patm=Patm*10**5
    
    rhogas = Patm/((8314.4598/molarmgas)*Tatm)
    Vbal= mpayload/(rhoatm-rhogas)*(buoyancyperc/100.)    
    mtot=(rhogas*Vbal+mpayload)/0.2    
    
    found=False
    while found==False:     

        Vbalnew=mtot/(rhoatm-rhogas)*(buoyancyperc/100.)
        mtotnew=(rhogas*Vbalnew+mpayload)/0.2

        if abs(mtotnew-mtot)<1e-3:
            found=True
        else:
            mtot=mtotnew
    mgas=rhogas*Vbalnew
    
    return(mtotnew, Vbalnew,mgas)

def fullbuoyancy(mgas, mtot, Vbal, expancruise):    
    import VenusAtmosphere as atm
        
    Vbalnew=Vbal*(100.-expancruise)/100.
    rhogasnew=mgas/Vbalnew

    #find what alt buoyancy is 100%
    #assumign the voume can increase by 10% but not decrease
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

import VenusAtmosphere as atm
hcruise = 50000 #km
Cl =1.5
mpayload = 90.
V=40.
hforce=hcruise/1000.
mtotnew, Vbalnew,mgas=ballooncalc(mpayload, hforce, 2.016, 100)
A=14

found=False
while not found:
    Tatm, Patm, rhoatm, GravAcc=atm.VenusAtmosphere30latitude(hforce)       
    Patm=Patm*10**5    
    rhogas = Patm/((8314.4598/2.016)*Tatm)
    Lift = GravAcc*Vbalnew*(rhoatm-rhogas)    
    altbuoy=fullbuoyancy(mgas,mtotnew,Vbalnew,10.)
    cord=(Vbalnew/2.4)**(1./3.)
    S=cord*A*cord
    LiftWing=1/2.*S*Cl*V**2*rhoatm
    if abs(Lift-LiftWing)<1e-1:
        found=True
    else:
        hforce-=1e-2
print hforce






