# -*- coding: utf-8 -*-
"""
Created on Wed May  4 16:57:58 2016

@author: MaxHenger
"""

import numpy as np

def isAscending(array):
    if len(array) == 0:
        raise ValueError("Expected at least one value in 'array'")
        I
    prev = array[0]
    
    for i in range(1, len(array)):
        if array[i] < prev:
            return False
            
        prev = array[i]
            
    return True
    
def isDescending(array):
    if len(array) == 0:
        raise ValueError("Expected at least one value in 'array'")
        
    prev = array[0]
    
    for i in range(1, len(array)):
        if array[i] > prev:
            return False
        
        prev = array[i]
        
    return True

def isArray(array):
    return isinstance(array, list) or isinstance(array, np.ndarray)
    
def DynViscosity(Temp,Viscinit=0.0000148,Tempinit=293.15,sutherland=240):
    """Returns the dynamic viscocity of CO2 as a function of temperature"""
    a = 0.555*Tempinit+sutherland
    b = 0.555*Temp + sutherland
    mu=Viscinit*(a/b)*np.power(Temp/Tempinit, 3/2)
    return mu
    
def scale_a():
    Rsp     = 192.5
    gamma   = 1.2941
    T=scale_height(0)[0]
    return np.sqrt(gamma*Rsp*T)
def scale_height(h):
    """Returns the atmospheric parameters as a function of altitude. Uses scale height method."""
    rho0    = 65.
    p0      = 92.1*10**5 
    Re      = 6052.*1000
    GM      = 0.32486*10**6*10**9
    H       = 15.9*1000
    Rsp     = 192.5
    
    rho = rho0*np.e**(-h/H)
    g = GM/(Re+h)**2
    p = p0*np.e**(-h/H)
    T = (p/(rho*Rsp)) # temperature is constant in scale height calculations
    return T,p,rho,g