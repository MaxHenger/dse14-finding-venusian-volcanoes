# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 11:51:01 2016

@author: MaxHenger
"""

import TrackCommon
import TrackIO

import numpy as np

class Settings:
    def __init__(self, version='iteration4'):
        version = version.lower()
        
        self.efficiencyCharging = 0.95 # Losses while charging the battery
        self.efficiencyDischarging = 0.95 # Losses while discharing the battery
        self.efficiencyPackingFactor = 0.85 # Packing factor of the solar cells
        self.efficiencyPower = 0.9 # Losses while consuming power
        self.efficiencyPropellers = 0.75 # Losses from shaft to thrust power

        self.PConsumption = 2342 # power consumption without propellers
        self.PPropeller = 32e3 # propeller power (NOT shaft power)
        self.reynoldsLength = 4.0 # value to use to determine the reynolds number
        
        # Venusian properties
        self.RVenus = 6051800
        self.omegaVenus = 2 * np.pi / (2802 * 24 * 60 * 60)
        self.muVenus = 0.32486e15
        self.fluxVenus = 2643
        
        # Track properties
        self.latitude = 0
        self.longitude = 0
        self.inclination = 0 # of the propellers, in radians
        self.speedOfSoundRatio = 0.6
        self.qInfMin = 200
        self.qInfMax = 1e10
        
        # Power properties
        self.specificWeightPanels = 0.84 # kg / m2
        self.batteryDepthOfDischarge = 80 # percent
        self.batterySafetyFactor = 0.25 # safety factor (ratio)
        
        self.WLander = 95

        if version == 'initial':
            # Aircraft properties
            self.W = 700 * 8.8
            self.S = 35
            self.SControl = 0

            # Lookup tables
            try:
                self.lookupCl, self.lookupCd = TrackIO.LoadAerodynamicData(
                    './data/aerodynamicPerformance/Cl.csv',
                    './data/aerodynamicPerformance/Cd.csv')
            except Exception as ex:
                self.lookupCl = None
                self.lookupCd = None
                print("Failed to load Cl/Cd:", ex)

            try:
                self.lowerBound, self.upperBound = TrackIO.LoadAscentGuides(
                    './optclimb_-60.0to20.0_0.0.dat', 4.0)
            except Exception as ex:
                print("Failed to load bounds:", ex)
                self.lowerBound = None
                self.upperBound = None

            self.alphaMin = -8.0
            self.alphaMax = 8.0

            
        elif version == 'iteration4':
            # Aircraft properties
            self.W = 858 * 8.8
            self.S = 40
            self.SControl = 14
            self.mBattery = 150
            
            # Lookup tables
            try:
                self.lookupCl, self.lookupCd = TrackIO.LoadAerodynamicReynoldsData(
                    './data/aerodynamicPerformance/v4/Cl.csv',
                    './data/aerodynamicPerformance/v4/Cd.csv')
            except Exception as ex:
                self.lookupCl = None
                self.lookupCd = None
                print("Failed to load Cl/Cd:", ex)
                
            self.lowerBound = None
            self.upperBound = None
            
            self.alphaMin = -8.0
            self.alphaMax = 8.0