# -*- coding: utf-8 -*-
"""
Created on Sun May 29 19:02:47 2016

@author: MaxHenger
"""

import numpy as np

import Atmosphere
import TrackDive
import TrackCommon

def OptimizeDive(heightUpper, heightTarget, vHorInitial, vVerInitial,
                 longitude, latitude, W, S, vHorTarget, dt,
                 lookupCl, lookupCd):
    # Retrieve ranges of angle of attack from the lookup tables
    alphaNew = lookupCl.getPoints()
    alphaNew = np.linspace(alphaNew[0][0], alphaNew[0][-1], 21)

    # Set initial values
    alpha = [0]
    vHor = [vHorInitial]
    vVer = [vVerInitial]
    height = [heightUpper]
    gamma = [0]

    # Set global variables
    atmosphere = Atmosphere.Atmosphere()

    # Start iterating
    while True:
        bias = 0.0

        while height[-1] > heightTarget:
            # Determine new flight variables
            vHorNew, vVerNew, gammaNew, hNew = TrackDive.Step(height[-1], alpha[-1],
                gamma[-1], vHor[-1], vVer[-1], longitude, latitude, W, S, alphaNew,
                dt, lookupCl, lookupCd, atmosphere)

            # Filter the valid solutions from the invalid ones
            iValid = []

            for i in range(0, len(alphaNew)):
                if gammaNew[i] > 0 and gammaNew[i] < np.pi / 2.0 and
