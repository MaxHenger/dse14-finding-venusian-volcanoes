# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 11:51:01 2016

@author: MaxHenger
"""

import TrackCommon

import numpy as np

class Settings:
    def __init__(self, version='initial'):
        version = version.lower()

        if version == 'initial':
            self.W = 700 * 8.8
            self.S = 35
            self.lookupCl, self.lookupCd = TrackCommon.LoadAerodynamicData(
                './data/aerodynamicPerformance/Cl.csv',
                './data/aerodynamicPerformance/Cd.csv')
            self.RVenus = 6051800
            self.omegaVenus = 2 * np.pi / (2802 * 24 * 60 * 60)
            self.muVenus = 0.32486e15
            self.latitude = 0
            self.longitude = 0
            self.inclination = 0 # of the propellers, in radians
