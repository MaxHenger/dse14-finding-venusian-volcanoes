"""
Created on Sun May 29 11:43:52 2016

@author: Mathijs

Thermal resistance functions 
"""

import numpy as np

def ThermalResistanceSphere(R_out,R_in,k):
    TR_sphere = (R_out - R_in) / (4*np.pi*k*R_in*R_out)
    return TR_sphere

def ThermalResistanceCylinder(R_out,R_in,k):
    TR_cylinder = np.log(R_out/R_in)/(2*np.pi*k)
    return TR_cylinder
    
def ThermalResistanceSlab(thickness,k,Area):
    TR_slab = thickness/(k*Area)
    return TR_slab