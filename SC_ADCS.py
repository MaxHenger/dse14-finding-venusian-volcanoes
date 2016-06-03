# -*- coding: utf-8 -*-
"""
Created on Mon May 23 14:08:34 2016

@author: Chaggai
"""
import numpy as np

#verified with formula from orbital mechanics
def disturbSRP(reflec, Ar_sc, dsun):
    #reflect - reflection coefficient of material
    # PSun is soalr pressure at 1 AU in N/m^2
    # Ar_sc is area of spacecraft visible to sun in m^2
    # dsun is distance between sun and satellite in m
    Psun = 4.5605*10**-6
    AU = 149598023*1000. #m
    F_SRP = (1.+reflec)*Psun*Ar_sc*(AU/dsun)**2.
    return(F_SRP)
    
def wheelsizing(Torquereq, maxomega):
    radius=0.001    
    found=False
    while found==False:        
        accel =0.051*radius**-2.021 #found with excel table relation (approximation)
        MoI = Torquereq/accel
        mass=(2.*MoI/(radius**2))
        accel2 = 16.817*mass**-0.927
        if abs(accel-accel2)<0.01:
            found=True
        else:
            radius+=0.0000001
        
    maxomega = np.deg2rad(maxomega*360)/60.
    maxh = MoI*maxomega
    return(mass, radius, maxh)
    
def momentdump(h,t,L):
    Ft = h/(t*L)
    return(Ft)
    
if __name__ == "__main__":
        
    Torquedist = 0.1 #N*m
    Marginfactor=1.5
    Torquereq = Torquedist*Marginfactor
    maxomega = 7500 #rpm
        
    
    print wheelsizing(Torquereq, maxomega)    
    

    
    


