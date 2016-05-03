# -*- coding: utf-8 -*-
"""
Created on Tue May 03 09:40:41 2016

@author: Julius
"""
import os
os.chdir("../../")  # change main directory to be able to import everything else


import solar
import VenusAtmosphere

# this is the concept design of a mxi between a blimp and an aircraft

Mass        = 2000  # [kg]
SurfaceArea = 20    # [m2]
thickness   = 1     # [m]
Volume      = thickness * SurfaceArea # [m3]

def func(boole):
    if boole:
        return boole
    return boole
