# -*- coding: utf-8 -*-
"""
Created on Tue May 10 19:46:30 2016

This file contains the constants that are used inside the Orbital modeling
code in Orbit.py (defining the Orbit class). This file does not need any 
editing if one is not attempting to edit or add a new model.

@author: Julius
"""

class OrbitVenus:
    def __init__(self):
        #Bulk Parameters
        self.Mass       = 4.8675 *10**24
        self.Volume     = 92.843 *10**10 *1000**3
        self.RadiusEq   = 6051.8 *1000
        self.RadiusMean = 6051.8 *1000
        self.RadiusPolar= 6051.8 *1000
        self.Flattening = 0.000
        self.Density    = 5243
        self.gravitationalAcceleration   = 8.87
        self.Mu         = 0.32486 *10**6 *10**9
        self.GM         = self.Mu
        self.Albedo     = 0.90
        #self.geometricAlbedo = 0.67
        #self.visualMagnitude = -4.4
        self.SolarIrradiance = 2601,3
        self.BlackBodyTemp = 184
        self.TopographicRange = 13 *1000
        self.J2         = 4.458 *10**-4
        
        #Orbit Parameters
        self.SemiMajor  = 108.21  *10**6
        self.SiderealOrbitPeriod = 224.701 *24*60*60
        self.TropicalOrbitPeriod = 224.695 *24*60*60
        self.Perihelion     = 107.48 *10**6
        self.Aperhelion     = 108.94 *10**6
        self.SynodicPeriod  = 583.92 *24*60*60
        self.VelocityMean   = 35.02 *1000
        self.VelocityMax    = 35.26 *1000
        self.VelocityMin    = 34.79 *1000
        self.Inclination    = 3.39 # deg
        self.Eccentricity   = 0.0067
        self.SiderealRotationalPeriod = -5832.8 *60*60
        self.DayLength      = 2802.0 *60*60
        self.Obliquity      = 177.36
        self.EquatorInclination = 2.64
        
        #North Pole of Rotation
        self.RightAscention = 272.76
        self.Declination = 67.16
        self.RefDate = ["12:00 UT 1 Jan 2000" , 2451545 ]
        
class OrbitEarth:
    def __init__(self):
        pass
    
    
class OrbitSun:
    def __init__(self):
        pass