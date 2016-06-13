# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 15:07:02 2016

This file is used to completely gauge the possible flight limits of our aircraft
using several predefined limits. Some experimental plotting is going on here as
well.

@author: MaxHenger
"""

import Atmosphere
import TrackCommon
import TrackAngleOfAttack
import TrackContour
import TrackSettings
import TrackClimbOptimize
import TrackPower
import TimeEstimator

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def GenerateCruiseMaps(axisHeight, axisDeltaV, latitude, longitude, relativeSeverity,
                       W, S, inclination, lookupCl, lookupCd, atm):
    # Generate the derivatives of Cl and Cd
    lookupdCldAlpha = lookupCl.getDerivative()
    lookupdCddAlpha = lookupCd.getDerivative()

    # Preallocate the resulting maps
    qInf = np.zeros([len(axisHeight), len(axisDeltaV)])
    vInf = np.zeros(qInf.shape)
    vInfOvera = np.zeros(qInf.shape)
    alpha = np.zeros(qInf.shape)
    thrust = np.zeros(qInf.shape)
    PReq = np.zeros(qInf.shape)

    velocityZonal = TrackCommon.AdjustSeverity(atm.velocityZonal(axisHeight, 
        latitude, longitude), relativeSeverity)
    density = TrackCommon.AdjustSeverity(atm.density(axisHeight, latitude, 
        longitude), relativeSeverity)
    speedOfSound = atm.speedOfSound(axisHeight, latitude, longitude)

    for iHeight in range(0, len(axisHeight)):
        # Calculate freesteram velocity and dynamic pressure
        vInfCur = abs(velocityZonal[iHeight] + axisDeltaV)
        vInf[iHeight, :] = vInfCur
        vInfOvera[iHeight, :] = vInfCur / speedOfSound[iHeight]
        qInf[iHeight, :] = 0.5 * density[iHeight] * np.power(vInfCur, 2.0)

        # Angle of attack and power required calculations are performed per
        # velocity value
        for iDeltaV in range(0, len(axisDeltaV)):
            alpha[iHeight, iDeltaV], thrust[iHeight, iDeltaV], _ = \
                TrackAngleOfAttack.AngleOfAttackThrustSteady(W, S,
                qInf[iHeight, iDeltaV], inclination, lookupCl, lookupdCldAlpha,
                lookupCd, lookupdCddAlpha)

            PReq[iHeight, iDeltaV] = thrust[iHeight, iDeltaV] * vInfCur[iDeltaV]

    return {"qInf": qInf, "vInf": vInf, "alpha": alpha,
            "thrust": thrust, "PReq": PReq, "vInfOvera": vInfOvera}
    

def IsEnoughPower(estimator, settings, atmosphere, heightMin, heightMax, 
                  latitude, longitude, severity, area, areaRatio, PReq,
                  lookupdCldAlpha, lookupdCddAlpha):
    for height in np.linspace(heightMin, heightMax, 100):
        # Figure out the power required to sustain cruise at the current height
        vZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(height, 
            latitude, longitude), severity)
        density = TrackCommon.AdjustSeverity(atmosphere.density(height,
            latitude, longitude), severity)
        
        aoa, thrust, valid = TrackAngleOfAttack.AngleOfAttackThrustSteady(
            settings.W, settings.S, 0.5 * density * vZonal**2.0,
            settings.inclination, settings.lookupCl, lookupdCldAlpha,
            settings.lookupCd, lookupdCddAlpha)
        
        if not valid:
            # cannot sustain cruise at this height
            continue
        
        currentPower = thrust * vZonal
        efficiencies = estimator.getPowerEfficiency(height, latitude, 
            longitude, aoa / 180.0 * np.pi, 0)
        
        if area * settings.fluxVenus * (efficiencies[0] + efficiencies[1] +
                                        efficiencies[2] * areaRatio) > currentPower:
            return True
    
    return False
    
def DetermineTracks(axisHeight, axisDeltaV, latitude, longitude, PReq,
                    axisSeverity):
    # For each of the provided severities generate cruise maps and figure out
    # where the aircraft can fly
    settings = TrackSettings.Settings()
    atmosphere = Atmosphere.Atmosphere()
    powerEstimator = TrackPower.TrackPower(settings, atmosphere)
    timeEstimator = TimeEstimator.TimeEstimator(len(axisSeverity))
    
    lookupdCldAlpha = settings.lookupCl.getDerivative()
    lookupdCddAlpha = settings.lookupCd.getDerivative()

    # Some bogus areas and area ratios as a lower bound for the heights at 
    # which sustained cruise is possible
    bogusArea = 40.0
    bogusAreaRatio = 2.5
    numUpdates = 10
    
    # Values that determine which points will be chosen as lower and upper
    # operational points
    numOperationalSeverities = 5 # In which sections to subdivide the operational severities
    operationalSeverityPercentage = 0.05 # percentage to remove from operational severity range
    numUpperOperationalPoints = 5
    lowerOperationalHeightBottom = 32e3
    lowerOperationalHeightTop = 50e3
    numLowerOperationalHeights = 4
    lowerDeltaVPerOperationalPoint = 10
    numPReqMinimisationPoints = 100
    
    severityBounds = [None, None]
    contourBounds = [None, None]

    print(TrackCommon.StringHeader("Finding valid atmospheric severities", 60))
    
    contour = None
    
    timeEstimator.startTiming()

    for iSeverity in range(0, len(axisSeverity)):
        timeEstimator.startIteration(iSeverity)
        
        maps = GenerateCruiseMaps(axisHeight, axisDeltaV, settings.latitude,
            settings.longitude, axisSeverity[iSeverity], settings.W,
            settings.S, settings.inclination, settings.lookupCl,
            settings.lookupCd, atmosphere)

        # Check if there are any valid cruise ranges that enable roughly
        # constant subsolar cruise (neglect rotation of Venus)
        contour = TrackContour.Contour()
        contour.combineData(axisDeltaV, axisHeight, [maps['qInf'], maps['vInf'],
            maps['alpha'], maps['PReq'], maps['vInfOvera']], [settings.qInfMin,
            5, settings.alphaMin, 0, 0], [settings.qInfMax, 1000,
            settings.alphaMax, PReq, settings.speedOfSoundRatio])

        rangesConstantSubsolar = contour.getVerticalRanges(0)

        isValid = False
        
        for iRange in range(0, len(rangesConstantSubsolar)):
            if IsEnoughPower(powerEstimator, settings, atmosphere,
                    rangesConstantSubsolar[iRange][0], rangesConstantSubsolar[iRange][1], 
                    settings.latitude, settings.longitude, axisSeverity[iSeverity],
                    bogusArea, bogusAreaRatio, PReq, lookupdCldAlpha, lookupdCddAlpha):
                isValid = True
                break
            
        if isValid:
            if severityBounds[0] == None:
                if iSeverity == 0:
                    print(' * Warning: First valid severity is at index 0. The ' +
                          'severity index can probably be further reduced.')
                          
                severityBounds[0] = axisSeverity[iSeverity]
                contourBounds[0] = contour
        elif severityBounds[0] != None and severityBounds[1] == None:
            severityBounds[1] = axisSeverity[iSeverity - 1]
            contourBounds[1] = contour
        
        timeEstimator.finishedIteration(iSeverity)
        
        print('.', end='')
        
        if (iSeverity + 1) % numUpdates == 0:
            print(' spent:', timeEstimator.getTotalElapsed(),
                  ', remaining:', timeEstimator.getEstimatedRemaining())
            
    if severityBounds[0] == None:
        raise ValueError("No valid severity value found")

    if severityBounds[1] == None:
        print(' * Warning: Final valid severity is at index', len(axisSeverity) - 1,
            ', the severity index can probably be further increased')
        severityBounds[1] = axisSeverity[-1]
        contourBounds[1] = contour
        print(' > Final valid severity:', round(severityBounds[1], 3))
    
    # Display some of the common minimum and maximum wind speeds at
    # various heights
#    print('')
#    names = ['minimum severity', 'maximum severity']
#
#    print(TrackCommon.StringHeader("Operational limits", 60))
#    
#    for iBound in range(0, len(severityBounds)):
#        print(" * Common windspeeds and densities at", names[iBound],
#              ":", severityBounds[iBound])
#        
#        for height in np.linspace(axisHeight[0] + 1, axisHeight[-1] - 1, 50):
#            speedBounds = contourBounds[iBound].getHorizontalRanges(height)
#            
#            if len(speedBounds) != 0:
#                minDeltaV = speedBounds[0][0]
#                maxDeltaV = speedBounds[-1][-1]
#                vZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(height, 
#                    latitude, longitude), severityBounds[iBound])
#                
#                print(TrackCommon.StringPad(" > h = ", height / 1e3, 1, 5) + 
#                      TrackCommon.StringPad(" km, vMin = ", minDeltaV + vZonal, 2, 6) + 
#                      TrackCommon.StringPad(" m/s, vMax = ", maxDeltaV + vZonal, 2, 6) + " m/s")
                
    print(TrackCommon.StringHeader("Determining operational cruise points", 60))
    severityDelta = operationalSeverityPercentage * (severityBounds[1] - severityBounds[0])
    operationalSeverities = np.linspace(severityBounds[0] + severityDelta,
        severityBounds[1] - severityDelta, numOperationalSeverities)
    
    finalData = []
    
    for iSeverity in range(0, len(operationalSeverities)):
        # Take the inclusive contours and see whether or not they cross the
        # 0 deltaV line. If they do seperate the contour into the negative-going
        # section and the positive going section.
#        print(' * Considering severity:', round(operationalSeverities[iSeverity], 3))
        maps = GenerateCruiseMaps(axisHeight, axisDeltaV, settings.latitude,
            settings.longitude, operationalSeverities[iSeverity], settings.W,
            settings.S, settings.inclination, settings.lookupCl, 
            settings.lookupCd, atmosphere)
        
        contour = TrackContour.Contour()
        contour.combineData(axisDeltaV, axisHeight, [maps['qInf'], maps['vInf'],
            maps['alpha'], maps['PReq'], maps['vInfOvera']], [settings.qInfMin,
            5, settings.alphaMin, 0, 0], [settings.qInfMax, 1000,
            settings.alphaMax, PReq, settings.speedOfSoundRatio])
        
        # Determine the upper operational points
        upperHeightMax = 0
        usableContour = None
        
        for iContour in range(0, contour.getNumContours()):
            curContour = contour.getContour(iContour)
            
            if len(curContour.getVerticalRanges(0)) == 0:
                continue
            
            if curContour.getMaxX() > upperHeightMax:
                upperHeightMax = curContour.getMaxX()
                usableContour = curContour
            
        upperHeightDelta = upperHeightMax / (2 * numUpperOperationalPoints)
        operationalDeltaV = np.linspace(upperHeightDelta, upperHeightMax -
            upperHeightDelta, numUpperOperationalPoints)
        
        # For each upper operational point find the minimum cruising power
        # needed.
        upperPReq = []
        upperHeight = []
        upperDeltaV = []

        for iOperationalDeltaV in range(0, numUpperOperationalPoints):
            #print(' > considering dV =', operationalDeltaV[iOperationalDeltaV])
            minimumPReq = PReq
            minimumPReqHeight = 0
            bestPowerRatio = 0
            curRange = usableContour.getVerticalRanges(operationalDeltaV[iOperationalDeltaV])
            
            #print(' > ranges:', curRange)
            for iRange in range(0, len(curRange)):
                curHeights = np.linspace(curRange[iRange][0], curRange[iRange][1],
                    numPReqMinimisationPoints)
                
                #print(' > heights:', curHeights)
                
                # Retrieve the local power required for each of the heights and
                # compare it to the current minimal value
                for iHeight in range(0, numPReqMinimisationPoints):
                    vZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(
                        curHeights[iHeight], settings.latitude, settings.longitude), 
                        operationalSeverities[iSeverity])
                    curSpeed = vZonal + operationalDeltaV[iOperationalDeltaV]
                    density = TrackCommon.AdjustSeverity(atmosphere.density(curHeights[iHeight],
                        settings.latitude, settings.longitude), 
                        operationalSeverities[iSeverity])
                    
                    aoa, thrust, valid = TrackAngleOfAttack.AngleOfAttackThrustSteady(
                        settings.W, settings.S, 0.5 * density * curSpeed**2.0,
                        settings.inclination, settings.lookupCl, lookupdCldAlpha,
                        settings.lookupCd, lookupdCddAlpha, tol=1e-5)
                    
                    if not valid:
                        continue
                    
                    curPower = thrust * curSpeed
                    
                    # Determine the power available per unit area
                    efficiencies = powerEstimator.getPowerEfficiency(curHeights[iHeight],
                        settings.latitude, settings.longitude, aoa / 180.0 * np.pi, 0)
                    
                    curRatio = (efficiencies[0] + efficiencies[1] + efficiencies[2] * 2.5) / curPower
                    
                    if curRatio > bestPowerRatio:
                        bestPowerRatio = curRatio
                        minimumPReq = curPower
                        minimumPReqHeight = curHeights[iHeight]

            upperPReq.append(minimumPReq)
            upperHeight.append(minimumPReqHeight)
            upperDeltaV.append(operationalDeltaV[iOperationalDeltaV])
            
#            print(TrackCommon.StringPad("at h = ", minimumPReqHeight / 1e3, 2, 6) + 
#                  TrackCommon.StringPad("km, dV = ", operationalDeltaV[iOperationalDeltaV], 2, 6) + 
#                  TrackCommon.StringPad("m/s, PReq = ", minimumPReq / 1e3, 3, 6) + " kW")
        
        operationalLowerHeight = np.linspace(lowerOperationalHeightBottom,
            lowerOperationalHeightTop, numLowerOperationalHeights)
        
        lowerDeltaV = []
        lowerHeight = []
        
        for iHeight in range(0, numLowerOperationalHeights):
            horizontalRanges = usableContour.getHorizontalRanges(operationalLowerHeight[iHeight])

            for iRange in range(0, len(horizontalRanges)):
                curRange = horizontalRanges[iRange]
                    
                if curRange[0] > 0:
                    continue
                
                if curRange[1] > 0:
                    curRange[1] = 0

                numPoints = int((curRange[1] - curRange[0]) / 
                    lowerDeltaVPerOperationalPoint) + 1
                
                if numPoints == 1:
                    lowerDeltaV.append((curRange[0] + curRange[1]) / 2.0)
                    lowerHeight.append(operationalLowerHeight[iHeight])
                else:
                    offset = (curRange[1] - curRange[0]) / (2.0 * numPoints)
                    
                    for i in range(0, numPoints):
                        lowerDeltaV.append(curRange[0] + offset + (curRange[1] - 
                            curRange[0] - 2 * offset) / (numPoints - 1) * i)
                        lowerHeight.append(operationalLowerHeight[iHeight])
            
#        # Print all results
#        for i in range(0, len(lowerDeltaV)):
#            print(TrackCommon.StringPad("at h = ", lowerHeight[i] / 1e3, 2, 6) +
#                  TrackCommon.StringPad("km, dV = ", lowerDeltaV[i], 2, 6) + " m/s")
#        
#        # Plot results
#        ax = plt.figure().add_subplot(111)
#        
#        for j in range(0, contour.getNumContours()):
#            curthing = contour.getContour(j)
#            includer = curthing.getIncluder()
#        
#            ax.plot(includer[:, 0], includer[:, 1])
#        ax.scatter(upperDeltaV, upperHeight, c=[1.0, 0.0, 0.0])
#        ax.scatter(lowerDeltaV, lowerHeight, c=[0.0, 1.0, 0.0])

        finalData.append([operationalSeverities[iSeverity],
                          upperDeltaV, upperHeight,
                          lowerDeltaV, lowerHeight])
        
    return finalData

#DetermineTracks(np.linspace(20, 80, 40) * 1e3, np.linspace(-50, 50, 40), 0, 0,
#                32e3, np.linspace(-2.0, 2.0, 50))