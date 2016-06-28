#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 11:36:39 2016

@author: MaxHenger
"""

import Atmosphere
import TrackCommon
import TrackIO
import TrackBounds
import TrackClimbOptimize
import TrackSettings
import TrackContour

import numpy as np
import matplotlib.pyplot as plt

def Generate(minHeight, maxHeight, numHeight, minSpeed, maxSpeed, numSpeed, 
             fileCl, fileCd, severity, resultFilename='mathijs.csv'):
    atm = Atmosphere.Atmosphere()
    lookupCl, lookupCd, = TrackIO.LoadAerodynamicReynoldsData(fileCl, fileCd)
    
    settings = TrackSettings.Settings()
    
    # Generate ascent and cruise maps
    axisHeight = np.linspace(minHeight, maxHeight, numHeight)
    axisDeltaV = np.linspace(minSpeed, maxSpeed, numSpeed)
    cruiseMaps = TrackBounds.GenerateCruiseMaps(axisHeight, axisDeltaV, 
        settings.latitude, settings.longitude, severity, settings.W, settings.S, 
        settings.inclination, lookupCl, lookupCd, atm, settings.reynoldsLength)
    
    ascentMaps = TrackClimbOptimize.GenerateAscentMaps(axisHeight, axisDeltaV,
        settings.W, settings.S, settings.inclination, lookupCl, lookupCd, atm,
        settings.qInfMin, settings.qInfMax, settings.alphaMin,
        settings.alphaMax, 0, 32e3, settings.latitude, settings.longitude,
        settings.reynoldsLength, severity)
    
    # Generate contours and find the one describing the mission space
    contour = TrackContour.Contour()
    contour.combineData(axisDeltaV, axisHeight, [cruiseMaps['qInf'], 
        cruiseMaps['alpha'], cruiseMaps['PReq'], cruiseMaps['vInfOvera'],
        ascentMaps['PReq'], ascentMaps['alpha'], ascentMaps['qInf']],
        [settings.qInfMin, settings.alphaMin, 0, -1, 0, settings.alphaMin,
         settings.qInfMin], [settings.qInfMax, settings.alphaMax, 32e3,
         settings.speedOfSoundRatio, 32e3, settings.alphaMax, settings.qInfMax])
    
    missionContour = None
    
    ax = plt.figure().add_subplot(111)
    for i in range(0, contour.getNumContours()):
        contourSection = contour.getContour(i)
        includer = contourSection.getIncluder()
        ax.plot(includer[:, 0], includer[:, 1])
        
        if len(contourSection.getVerticalRanges(0)) != 0:
            missionContour = contourSection
            break
        
    if missionContour is None:
        raise ValueError("Failed to find mission bounding contour")
    
    # Prepare arrays for result storage
    density = np.zeros(axisHeight.shape)
    vZonal = np.zeros(axisHeight.shape)
    minVelocity = np.zeros(axisHeight.shape)
    maxVelocity = np.zeros(axisHeight.shape)
    stallSpeedLander = np.zeros(axisHeight.shape)
    stallSpeedLanderless = np.zeros(axisHeight.shape)
    thrustMinVelocity = np.zeros(axisHeight.shape)
    thrustMaxVelocity = np.zeros(axisHeight.shape)
    gravityAcceleration = np.zeros(axisHeight.shape)
    
    numTestSpeed = 32
    
    for iHeight in range(len(axisHeight)):
        density[iHeight] = TrackCommon.AdjustSeverity(atm.density(
            axisHeight[iHeight], settings.latitude, settings.longitude), severity)
        vZonal[iHeight] = TrackCommon.AdjustSeverity(atm.velocityZonal(
            axisHeight[iHeight], settings.latitude, settings.longitude), severity)
        gravityAcceleration[iHeight] = atm.gravitationalAcceleration(axisHeight[iHeight])
        
        ranges = missionContour.getHorizontalRanges(axisHeight[iHeight])
        
        if len(ranges) == 0:
            minVelocity[iHeight] = -1.0
            maxVelocity[iHeight] = -1.0
            thrustMinVelocity[iHeight] = -1.0
            thrustMaxVelocity[iHeight] = -1.0
        else:
            minVelocity[iHeight] = ranges[0][0] + vZonal[iHeight]
            maxVelocity[iHeight] = ranges[-1][1] + vZonal[iHeight]
            thrustMinVelocity[iHeight] = 32e3 / minVelocity[iHeight]
            thrustMaxVelocity[iHeight] = 32e3 / maxVelocity[iHeight]
            
        # Iterate towards the stall speed
        weights = [settings.W, settings.W - settings.WLander]

        for iWeight in range(0, 2):
            vLower = 0
            vUpper = vZonal[iHeight] + maxVelocity[iHeight]
            
            for iIt in range(0, numTestSpeed):
                vCenter = (vLower + vUpper) / 2
                Re = TrackCommon.AdjustSeverity(atm.reynoldsNumber(axisHeight[iHeight],
                    settings.latitude, settings.longitude, vCenter,
                    settings.reynoldsLength), severity)
                
                curCl = lookupCl.find(settings.alphaMax, Re)
                vStall = np.sqrt(2 * weights[iWeight] / (density[iHeight] * 
                    settings.S * curCl))
                
                if vStall > vCenter:
                    # Stall speed is larger then the assumed speed for the
                    # therefrom deduced reynolds number: increase assumed speed
                    vLower = vCenter
                else:
                    vUpper = vCenter
                
            vCenter = (vLower + vUpper) / 2

            if iWeight == 0:
                stallSpeedLander[iHeight] = vCenter
            else:
                stallSpeedLanderless[iHeight] = vCenter

    # Store all data in a file
    np.savetxt(resultFilename, np.transpose([axisHeight, density, vZonal,
        minVelocity, maxVelocity, stallSpeedLander, stallSpeedLanderless,
        thrustMinVelocity, thrustMaxVelocity, gravityAcceleration]),
        '%8.8f', ';', header='height [m]; density [kg/m3]; zonal wind [m/s]; ' +
        'climb velocity (min) [m/s]; cruise velocity (max) [m/s]; ' +
        'stall speed [m/s]; stall speed without lander [m/s]; ' +
        'thrust at min velocity [N]; thrust at max velocity [N]; ' +
        'gravitational acceleration [m/s2]')
                    
Generate(32e3, 70e3, 39, -80, 80, 95, './data/aerodynamicPerformance/v4/Cl.csv',
    './data/aerodynamicPerformance/v4/Cd.csv', -1.0, 'mathijs_-1.0.csv')
Generate(32e3, 70e3, 39, -80, 80, 95, './data/aerodynamicPerformance/v4/Cl.csv',
    './data/aerodynamicPerformance/v4/Cd.csv', 0.0, 'mathijs_0.0.csv')
Generate(32e3, 70e3, 39, -80, 80, 95, './data/aerodynamicPerformance/v4/Cl.csv',
    './data/aerodynamicPerformance/v4/Cd.csv', 1.0, 'mathijs_1.0.csv')