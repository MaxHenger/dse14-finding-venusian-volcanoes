"""
Created on Sun May 29 11:43:52 2016

@author: Mathijs

Thermal resistance functions 
"""

import numpy as np

'''
Inputs:
R_out = outer radius (m)
R_in = inner radius (m)
k = conductivity of material [W/(m*K)]
'''

def ThermalResistanceSphere(R_out,R_in,k):
    TR_sphere = (R_out - R_in) / (4*np.pi*k*R_in*R_out)
    return TR_sphere
    
'''
Inputs:
R_out = outer radius (m)
R_in = inner radius (m)
k = conductivity of material [W/(m*K)]
L = length of cylinder (m)
'''

def ThermalResistanceCylinder(R_out,R_in,k,L):
    TR_cylinder = np.log(R_out/R_in)/(2*np.pi*k*L)
    return TR_cylinder
    
'''
Inputs:
thickness = thickness of material [m]
k = conductivity of material [W/(m*K)]
Area = area of conduction [m^2]    
'''

def ThermalResistanceSlab(thickness,k,Area):
    TR_slab = thickness/(k*Area)
    return TR_slab