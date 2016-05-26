# -*- coding: utf-8 -*-
"""
Created on Mon May 23 14:08:34 2016

@author: Chaggai
"""
import numpy as np

#verified with formula from orbital mechanics
def disturbSRP(reflec,Psun, Ar_sc, dsun):
    #reflect - reflection coefficient of material
    # PSun is soalr pressure at 1 AU in N/m^2
    # Ar_sc is area of spacecraft visible to sun in m^2
    # dsun is distance between sun and satellite in m
    AU = 149598023*1000. #m
    F_SRP = (1.+reflec)*Psun*Ar_sc*(AU/dsun)**2.
    return(F_SRP)
    

P_sun1AU = 4.5605*10**-6 #N/m^2


