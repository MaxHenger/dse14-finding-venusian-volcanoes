# -*- coding: utf-8 -*-
"""
Created on Tue May 10 19:46:30 2016

This file contains the constants that are used inside the gravity modeling
code in Gravity.py (defining the Gravity class). This file does not need any 
editing if one is not attempting to edit or add a new model.

I know this file is very ... sparse, however if we want to develop this further
it will serve as a good basis.

@author: Julius
"""
import numpy as np

class GravityVenus:
    def __init__(self):
        self.Mu         = 0.32486 *10**6 *10**9
        self.RadiusMean = 6051.8 *1000
        self.J          = np.zeros((3,3))
        self.lam        = np.zeros((3,3))
        
        self.J[2][0]    = 4.458 *10**-4  # the J2 effect
        
class GravityEarth:
    def __init__(self):
        self.J          = np.zeros((7,7))
        self.lam        = np.zeros((7,7))
        self.J[2][0]    = 1082.6267 *10**-6 # the J2 effect
        self.J[3][0]    = -2.5327 *10**-6
        self.J[4][0]    = -1.6196 *10**-6 
        self.J[5][0]    = -0.2273 *10**-6 
    
class GravitySun:
    def __init__(self):
        pass