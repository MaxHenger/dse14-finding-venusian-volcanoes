# -*- coding: utf-8 -*-
"""
Created on Tue May 10 19:08:39 2016

This file contains the definition of an orbit as a class with methods to
calculate the position, solar flux, ect. of the sattelite.



In general all functions accept and return values in SI units. An exception is
angles, these are entered in degrees. Also attempt to call the functions with
arrays as much as possible to reduce unneccesary calculations.

The general use of this file is as following:
    
    import Orbit
    orb = Orbit.Orbit('preliminary')
    orb.create(SemiMajor,Eccentricity,Inclination,AscentionAngle)
    orb.place(LocationAngle,Time)
    orb.where(Time)

    

@author: Julius
"""

import matplotlib.pyplot as plt
import numpy as np
import Utility as util
import Atmosphere
import OrbitConstants as orb_const

class Orbit:
    """Orbit class defining the orbit of the spacecraft"""
    def __init__(self,grav,SemiMajor,Eccentricity,Inclination,AscentionAngle,ArgumentPeri,TrueAnomaly,Time=0,Orbiting="Venus"):
        if Orbiting=="Venus":
            self.constants=orb_const.OrbitVenus()
        elif Orbiting=="Earth":
            self.constants=orb_const.OrbitEarth()
        elif Orbiting=="Sun":
            self.constants=orb_const.OrbitSun()
        else:
            raise ValueError("Orbiting object not found.")
        
        self.init={}
        self.parameters={}
        self.grav=grav
        self.__create__(SemiMajor,Eccentricity,Inclination,AscentionAngle,ArgumentPeri)
        self.__place__(TrueAnomaly,Time)
        
    def __create__(self,SemiMajor,Eccentricity,Inclination,RightAscention,ArgumentPeri):
        self.SemiMajor      = SemiMajor
        self.Eccentricity   = Eccentricity
        assert Inclination>=0 and Inclination<180
        self.Inclination    = Inclination
        assert RightAscention>=0 and RightAscention<360
        self.RightAscention = RightAscention
        assert ArgumentPeri>=0 and ArgumentPeri<360
        self.ArgumentPeri   = ArgumentPeri
        
        self.Period         = 2*np.pi*( self.SemiMajor**3 / self.constants.Mu  )**0.5 
        self.Periapsis      = self.SemiMajor*(1-self.Eccentricity**2)/(1+self.Eccentricity*np.cos(np.deg2rad(0)))
        if self.Eccentricity<1:
            self.Apoapsis   = self.SemiMajor*(1-self.Eccentricity**2)/(1+self.Eccentricity*np.cos(np.deg2rad(180)))
        else:
            self.Apoapsis   = nan
        self.Energy         = - self.constants.Mu/(2*self.SemiMajor)
        
    def __place__(self,TrueAnomaly,Time):
        """Location is given in the angle from perapsis. Time is the baseline where all the other time is relative to"""
        self.init["LocationAngle"]=TrueAnomaly
        self.init["Time"]=Time
        
    def status(self,Time):
        pass
        
    def show(self,Time):
        pass
        
    def __where__(self,Time):
        """Returns location in orbital and spherical coordinates """
        r = self.parametes
        
    def __velocity__(self,Time):
        """Returns velocity in direction of orbit"""
        
    #def simulate(self,Tbegin=0,Tend,Tstep): # generates an array of Time, spehrical coords, velocity,  
        
        