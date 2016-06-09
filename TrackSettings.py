# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 11:51:01 2016

@author: MaxHenger
"""

import TrackCommon
import TrackIO

import numpy as np

class Settings:
    def __init__(self, version='initial'):
        version = version.lower()

        if version == 'initial':
            # Aircraft properties
            self.W = 700 * 8.8
            self.S = 35

            self.efficiencyCharging = 0.98 # Losses while charging the battery
            self.efficiencyDischarging = 0.98 # Losses while discharing the battery
            self.efficiencyPackingFactor = 0.85 # Packing factor of the solar cells
            self.efficiencyPower = 0.9 # Losses while consuming power
            self.efficiencyPropellers = 0.8 # Losses from shaft to thrust power

            self.PConsumption = 1210 # power consumption without propellers
            # 1 kW comes from the control surfaces

            # Lookup tables
            self.lookupCl, self.lookupCd = TrackIO.LoadAerodynamicData(
                './data/aerodynamicPerformance/Cl.csv',
                './data/aerodynamicPerformance/Cd.csv')
            self.lowerBound, self.upperBound = TrackIO.LoadAscentGuides(
                './optclimb_-60.0to20.0_0.0.dat')

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
            self.alphaMin = -8.0
            self.alphaMax = 8.0

            # Power properties
            self.specificWeightPanels = 0.84 # kg / m2
            self.batteryDepthOfDischarge = 80 # percent
            self.batterySafetyFactor = 0.25 # safety factor (ratio)
