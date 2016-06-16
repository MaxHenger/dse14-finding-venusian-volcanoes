af# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 00:58:02 2016

@author: MaxHenger
"""

import Atmosphere
import TrackCommon
import TrackDiveOptimize
import TrackAcceleratingOptimize
import TrackClimbOptimize
import TrackSettings
import TrackStorage
import TrackPower
import TimeEstimator
import AircraftBatteries

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as scp_int

def GroundRotation(time, height, vVer, settings):
    y = vVer / (settings.RVenus + height) - settings.omegaVenus
    return scp_int.trapz(y, time)

def SetAxisColors(ax, color):
    for tick in ax.get_yticklabels():
        tick.set_color(color)

def StitchTracks(preDiveHeight, preDiveVHor, postDiveHeight, postDiveVHor,
                 postDiveLoiter, preAscentVHor, postAscentVHor,
                 PReqMin, PReqMax, dt, settings, severity, plotResult=True,
                 saveResult=True):
    # For less verbose typing
    W = settings.W
    S = settings.S
    latitude = settings.latitude
    longitude = settings.longitude
    inclination = settings.inclination

    # Load all required data
    atm = Atmosphere.Atmosphere()

    # Start by performing a dive
    timeDive, heightDive, vHorDive, vVerDive, vInfDive, alphaDive, gammaDive = \
        TrackDiveOptimize.OptimizeDive(preDiveHeight, postDiveHeight,
        preDiveVHor, 0, longitude, latitude, W, S, postDiveVHor, 0,
        settings.speedOfSoundRatio, dt, settings.lookupCl, settings.lookupCd,
        severity, plotResults=False, storeResults=False)
    powerDive = np.zeros([len(timeDive)])

    timeEndDive = timeDive[-1] + dt

    # Change the speed to the desired one, use the final values from the dive
    # as input
    vZonalAcc1 = TrackCommon.AdjustSeverity(atm.velocityZonal(heightDive[-1],
        latitude, longitude), severity)

    timeAcc1, vHorAcc1, alphaAcc1, powerAcc1, vHorAvgAcc1, powerAvgAcc1 = \
        TrackAcceleratingOptimize.OptimizeAccelerating(heightDive[-1],
        vHorDive[-1], 0, alphaDive[-1], longitude, latitude, W, S, postDiveVHor,
        inclination, dt, PReqMin, PReqMax, settings.speedOfSoundRatio,
        settings.lookupCl, settings.lookupCd, severity, plotResults=False,
        storeResults=False)
    heightAcc1 = np.repeat(heightDive[-1], len(timeAcc1))
    vVerAcc1 = np.zeros([len(timeAcc1)])
    vInfAcc1 = vZonalAcc1 + vHorAcc1
    gammaAcc1 = np.zeros([len(timeAcc1)])

    print('acc1: ' + TrackCommon.StringPad('vZonal = ', vZonalAcc1, 2, 5) +
          TrackCommon.StringPad(' m/s, vHorAcc1[0] = ', vHorAcc1[0], 2, 5) +
          TrackCommon.StringPad(' m/s, vHorAcc1[-1] = ', vHorAcc1[-1], 2, 5) +
          TrackCommon.StringPad(' m/s, vInfAcc1[0] = ', vInfAcc1[0], 2, 5) +
          TrackCommon.StringPad(' m/s, vInfAcc1[-1] = ', vInfAcc1[-1], 2, 5) +
          TrackCommon.StringPad(' m/s, vHorAvgAcc1 = ', vHorAvgAcc1, 2, 5) + ' m/s')

    timeEndAcc1 = timeAcc1[-1] + timeEndDive + dt

    # Performing loiter. Appending these arrays (inefficiently) because one day
    # I might actually implement some proper loitering
    postDiveLoiterNum = int(postDiveLoiter / dt)
    if postDiveLoiterNum == 0:
        # waste at least one dt to make writing code easier
        postDiveLoiterNum = 1

    timeLoiter1 = np.linspace(0, (postDiveLoiterNum - 1) * dt, postDiveLoiterNum)
    heightLoiter1 = np.repeat(heightAcc1[-1], postDiveLoiterNum)
    vHorLoiter1 = np.repeat(vHorAvgAcc1, postDiveLoiterNum)
    vVerLoiter1 = np.zeros([postDiveLoiterNum])
    vInfLoiter1 = np.repeat(vInfAcc1[-1], postDiveLoiterNum)
    alphaLoiter1 = np.repeat(alphaAcc1[-1], postDiveLoiterNum)
    gammaLoiter1 = np.repeat(gammaAcc1[-1], postDiveLoiterNum)
    powerLoiter1 = np.repeat(powerAvgAcc1, postDiveLoiterNum)

    timeEndLoiter1 = timeLoiter1[-1] + timeEndAcc1 + dt

    # Perform the second speed change before initiating the climb
    vZonalAcc2 = TrackCommon.AdjustSeverity(atm.velocityZonal(heightLoiter1[-1],
        latitude, longitude), severity)

    timeAcc2, vHorAcc2, alphaAcc2, powerAcc2, vHorAvgAcc2, powerAvgAcc2 = \
        TrackAcceleratingOptimize.OptimizeAccelerating(heightLoiter1[-1],
        vHorLoiter1[-1], powerLoiter1[-1], alphaLoiter1[-1], longitude,
        latitude, W, S, preAscentVHor, inclination, dt, PReqMin, PReqMax,
        settings.speedOfSoundRatio, settings.lookupCl, settings.lookupCd,
        severity, plotResults=False, storeResults=False)
    heightAcc2 = np.repeat(heightLoiter1[-1], len(timeAcc2))
    vVerAcc2 = np.zeros([len(timeAcc2)])
    vInfAcc2 = vZonalAcc2 + vHorAcc2
    gammaAcc2 = np.zeros([len(timeAcc2)])

    print('acc2: ' + TrackCommon.StringPad('vZonal = ', vZonalAcc2, 2, 5) +
          TrackCommon.StringPad(' m/s, vHorAcc2[0] = ', vHorAcc2[0], 2, 5) +
          TrackCommon.StringPad(' m/s, vHorAcc2[-1] = ', vHorAcc2[-1], 2, 5) +
          TrackCommon.StringPad(' m/s, vInfAcc2[0] = ', vInfAcc2[0], 2, 5) +
          TrackCommon.StringPad(' m/s, vInfAcc2[-1] = ', vInfAcc2[-1], 2, 5) +
          TrackCommon.StringPad(' m/s, vHorAvgAcc2 = ', vHorAvgAcc2, 2, 5) + ' m/s')

    timeEndAcc2 = timeAcc2[-1] + timeEndLoiter1 + dt

    # Start the ascent
    heightClimbQuit = heightAcc1[-1] - 0.1 * (preDiveHeight - heightAcc1[-1])

    timeClimb, heightClimb, vHorClimb, vVerClimb, vInfClimb, alphaClimb, gammaClimb = \
        TrackClimbOptimize.OptimizeClimb(heightAcc2[-1], preDiveHeight,
        heightClimbQuit, vHorAvgAcc2, vVerAcc2[-1], longitude, latitude, W, S,
        postAscentVHor, 0, PReqMax, settings.speedOfSoundRatio, inclination, dt,
        settings.lookupCl, settings.lookupCd, severity,
        lookupBoundLowerVInf=settings.lowerBound,
        lookupBoundUpperVInf=settings.upperBound, plotResults=False,
        storeResults=False)
    powerClimb = np.repeat(PReqMax, len(timeClimb))

    timeEndClimb = timeClimb[-1] + timeEndAcc2 + dt

    # Accelerate to post ascent horizontal velocity
    vZonalAcc3 = TrackCommon.AdjustSeverity(atm.velocityZonal(heightClimb[-1],
        latitude, longitude), severity)

    timeAcc3, vHorAcc3, alphaAcc3, powerAcc3, vHorAvgAcc3, powerAvgAcc3 = \
        TrackAcceleratingOptimize.OptimizeAccelerating(heightClimb[-1],
        vHorClimb[-1], powerClimb[-1], alphaClimb[-1], longitude, latitude,
        W, S, postAscentVHor, inclination, dt, PReqMin, PReqMax,
        settings.speedOfSoundRatio, settings.lookupCl, settings.lookupCd,
        severity, plotResults=False, storeResults=False)
    heightAcc3 = np.repeat(heightClimb[-1], len(timeAcc3))
    vVerAcc3 = np.zeros([len(timeAcc3)])
    vInfAcc3 = vZonalAcc3 + vHorAcc3
    gammaAcc3 = np.zeros([len(timeAcc3)])

    print('acc3: ' + TrackCommon.StringPad('vZonal = ', vZonalAcc3, 2, 5) +
          TrackCommon.StringPad(' m/s, vHorAcc3[0] = ', vHorAcc3[0], 2, 5) +
          TrackCommon.StringPad(' m/s, vHorAcc3[-1] = ', vHorAcc3[-1], 2, 5) +
          TrackCommon.StringPad(' m/s, vInfAcc3[0] = ', vInfAcc3[0], 2, 5) +
          TrackCommon.StringPad(' m/s, vInfAcc3[-1] = ', vInfAcc3[-1], 2, 5) +
          TrackCommon.StringPad(' m/s, vHorAvgAcc3 = ', vHorAvgAcc3, 2, 5) + ' m/s')

    timeEndAcc3 = timeAcc3[-1] + timeEndClimb + dt

    # Already perform the post-loiter acceleration before diving. This is used
    # to estimate the time needed at loitering to end up at the same subsolar
    # point in the end
    vZonalAcc4 = TrackCommon.AdjustSeverity(atm.velocityZonal(heightAcc3[-1],
        latitude, longitude), severity)

    timeAcc4, vHorAcc4, alphaAcc4, powerAcc4, vHorAvgAcc4, powerAvgAcc4 = \
        TrackAcceleratingOptimize.OptimizeAccelerating(heightAcc3[-1],
        vHorAvgAcc3, powerAvgAcc3, alphaAcc3[-1], longitude, latitude,
        W, S, preDiveVHor, inclination, dt, PReqMin, PReqMax,
        settings.speedOfSoundRatio, settings.lookupCl, settings.lookupCd,
        severity, plotResults=False, storeResults=False)
    heightAcc4 = np.repeat(heightAcc3[-1], len(timeAcc4))
    vVerAcc4 = np.zeros([len(timeAcc4)])
    vInfAcc4 = vZonalAcc4 + vHorAcc4
    gammaAcc4 = np.zeros([len(timeAcc4)])

    print('acc4: ' + TrackCommon.StringPad('vZonal = ', vZonalAcc4, 2, 5) +
          TrackCommon.StringPad(' m/s, vHorAcc4[0] = ', vHorAcc4[0], 2, 5) +
          TrackCommon.StringPad(' m/s, vHorAcc4[-1] = ', vHorAcc4[-1], 2, 5) +
          TrackCommon.StringPad(' m/s, vInfAcc4[0] = ', vInfAcc4[0], 2, 5) +
          TrackCommon.StringPad(' m/s, vInfAcc4[-1] = ', vInfAcc4[-1], 2, 5) +
          TrackCommon.StringPad(' m/s, vHorAvgAcc4 = ', vHorAvgAcc4, 2, 5) + ' m/s')

    # Determine the upper-height loiter time required to reach the same subsolar
    # point again.
    rotationCovered = 0.0

    for t, h, v in [(timeDive, heightDive, vHorDive),
                    (timeAcc1, heightAcc1, vHorAcc1),
                    (timeLoiter1, heightLoiter1, vHorLoiter1),
                    (timeAcc2, heightAcc2, vHorAcc2),
                    (timeClimb, heightClimb, vHorClimb),
                    (timeAcc3, heightAcc3, vHorAcc3),
                    (timeAcc4, heightAcc4, vHorAcc4)]:
        rotationCovered += GroundRotation(np.asarray(t), np.asarray(h),
                                          np.asarray(v), settings)

    timeLoiterUp = ((settings.RVenus + heightAcc3[-1]) / (settings.omegaVenus *
        (settings.RVenus + heightAcc3[-1]) - vHorAvgAcc3) * rotationCovered)

    if timeLoiterUp < 0:
        raise ValueError("Upper loiter time is smaller than 0:", timeLoiterUp)

    postClimbLoiterNum = int(timeLoiterUp / dt)

    print(' > spending', round(timeLoiterUp, 1), 's at top (', postClimbLoiterNum, 'iterations )')

    if postClimbLoiterNum == 0:
        postClimbLoiterNum = 1

    timeLoiter2 = np.linspace(0, (postClimbLoiterNum - 1) * dt, postClimbLoiterNum)
    heightLoiter2 = np.repeat(heightAcc3[-1], postClimbLoiterNum)
    vHorLoiter2 = np.repeat(vHorAvgAcc3, postClimbLoiterNum)
    vVerLoiter2 = np.zeros([postClimbLoiterNum])
    vInfLoiter2 = np.repeat(vInfAcc3[-1], postClimbLoiterNum)
    alphaLoiter2 = np.repeat(alphaAcc3[-1], postClimbLoiterNum)
    gammaLoiter2 = np.repeat(gammaAcc3[-1], postClimbLoiterNum)
    powerLoiter2 = np.repeat(powerAvgAcc3, postClimbLoiterNum)

    timeEndLoiter2 = timeLoiter2[-1] + timeEndAcc3 + dt
    timeEndAcc4 = timeAcc4[-1] + timeEndLoiter2 + dt

    # Compound all data
    # - create empty arrays
    timeTotal = []
    heightTotal = []
    vHorTotal = []
    vVerTotal = []
    vInfTotal = []
    alphaTotal = []
    gammaTotal = []
    powerTotal = []

    # - glue all data together in a plottable manner
    timeTotal.extend(timeDive)
    heightTotal.extend(heightDive)
    vHorTotal.extend(vHorDive)
    vVerTotal.extend(vVerDive)
    vInfTotal.extend(vInfDive)
    alphaTotal.extend(alphaDive)
    gammaTotal.extend(gammaDive)
    powerTotal.extend(powerDive)

    timeTotal.extend(np.asarray(timeAcc1) + timeEndDive)
    heightTotal.extend(heightAcc1)
    vHorTotal.extend(vHorAcc1)
    vVerTotal.extend(vVerAcc1)
    vInfTotal.extend(vInfAcc1)
    alphaTotal.extend(alphaAcc1)
    gammaTotal.extend(gammaAcc1)
    powerTotal.extend(powerAcc1)

    timeTotal.extend(np.asarray(timeLoiter1) + timeEndAcc1)
    heightTotal.extend(heightLoiter1)
    vHorTotal.extend(vHorLoiter1)
    vVerTotal.extend(vVerLoiter1)
    vInfTotal.extend(vInfLoiter1)
    alphaTotal.extend(alphaLoiter1)
    gammaTotal.extend(gammaLoiter1)
    powerTotal.extend(powerLoiter1)

    timeTotal.extend(np.asarray(timeAcc2) + timeEndLoiter1)
    heightTotal.extend(heightAcc2)
    vHorTotal.extend(vHorAcc2)
    vVerTotal.extend(vVerAcc2)
    vInfTotal.extend(vInfAcc2)
    alphaTotal.extend(alphaAcc2)
    gammaTotal.extend(gammaAcc2)
    powerTotal.extend(powerAcc2)

    timeTotal.extend(np.asarray(timeClimb) + timeEndAcc2)
    heightTotal.extend(heightClimb)
    vHorTotal.extend(vHorClimb)
    vVerTotal.extend(vVerClimb)
    vInfTotal.extend(vInfClimb)
    alphaTotal.extend(alphaClimb)
    gammaTotal.extend(gammaClimb)
    powerTotal.extend(powerClimb)

    timeTotal.extend(np.asarray(timeAcc3) + timeEndClimb)
    heightTotal.extend(heightAcc3)
    vHorTotal.extend(vHorAcc3)
    vVerTotal.extend(vVerAcc3)
    vInfTotal.extend(vInfAcc3)
    alphaTotal.extend(alphaAcc3)
    gammaTotal.extend(gammaAcc3)
    powerTotal.extend(powerAcc3)

    timeTotal.extend(np.asarray(timeLoiter2) + timeEndAcc3)
    heightTotal.extend(heightLoiter2)
    vHorTotal.extend(vHorLoiter2)
    vVerTotal.extend(vVerLoiter2)
    vInfTotal.extend(vInfLoiter2)
    alphaTotal.extend(alphaLoiter2)
    gammaTotal.extend(gammaLoiter2)
    powerTotal.extend(powerLoiter2)

    timeTotal.extend(np.asarray(timeAcc4) + timeEndLoiter2)
    heightTotal.extend(heightAcc4)
    vHorTotal.extend(vHorAcc4)
    vVerTotal.extend(vVerAcc4)
    vInfTotal.extend(vInfAcc4)
    alphaTotal.extend(alphaAcc4)
    gammaTotal.extend(gammaAcc4)
    powerTotal.extend(powerAcc4)

    # - convert to numpy arrays for data manipulation
    timeTotal = np.asarray(timeTotal)
    heightTotal = np.asarray(heightTotal)
    vHorTotal = np.asarray(vHorTotal)
    vVerTotal = np.asarray(vVerTotal)
    vInfTotal = np.asarray(vInfTotal)
    alphaTotal = np.asarray(alphaTotal)
    gammaTotal = np.asarray(gammaTotal)
    powerTotal = np.asarray(powerTotal)

    # Start plotting
    if plotResult:
        fig = plt.figure()
        axHeight = fig.add_subplot(311)
        axVelCom = fig.add_subplot(312, sharex=axHeight)
        axVelInf = axVelCom.twinx()
        axAlpha = fig.add_subplot(313, sharex=axHeight)
        axGamma = axAlpha.twinx()
    
        axHeight.plot(timeTotal, heightTotal / 1e3, 'g')
        axHeight.set_xlabel('time [s]')
        axHeight.set_ylabel('height [km]')
        axHeight.grid(True)
    
        axVelCom.plot(timeTotal, vHorTotal, 'r', label='vHor')
        axVelCom.plot(timeTotal, vVerTotal, 'r--', label='vVer')
        axVelInf.plot(timeTotal, vInfTotal, 'g', label='vInf')
    
        SetAxisColors(axVelCom, 'r')
        SetAxisColors(axVelInf, 'g')
    
        axVelCom.set_xlabel('time [s]')
        axVelCom.set_ylabel('velocity [m/s]')
        axVelInf.set_ylabel('velocity [m/s]')
        axVelCom.legend()
        axVelCom.grid(True)
    
        axAlpha.plot(timeTotal, alphaTotal, 'r', label='alpha')
        axGamma.plot(timeTotal, gammaTotal, 'g', label='gamma')
    
        SetAxisColors(axAlpha, 'r')
        SetAxisColors(axGamma, 'g')
    
        axAlpha.set_xlabel('time [s]')
        axAlpha.set_ylabel('alpha [deg]')
        axGamma.set_ylabel('gamma [deg]')
    
        fig.suptitle('Flight parameters')

    # Calculate the covered ground
    groundCovered = TrackCommon.CumulativeSimps(vHorTotal * settings.RVenus /
        (settings.RVenus + heightTotal), timeTotal)
    solarGroundCovered = groundCovered - settings.omegaVenus * settings.RVenus * timeTotal

    # Store all data
    finalFilename = None
    
    if saveResult:
        file = TrackStorage.DataStorage()
        file.addVariable('timeTotal', timeTotal)
        file.addVariable('heightTotal', heightTotal)
        file.addVariable('vHorTotal', vHorTotal)
        file.addVariable('vVerTotal', vVerTotal)
        file.addVariable('vInfTotal', vInfTotal)
        file.addVariable('alphaTotal', alphaTotal)
        file.addVariable('gammaTotal', gammaTotal)
        file.addVariable('powerTotal', powerTotal)
        file.addVariable('latitude', latitude)
        file.addVariable('longitude', longitude)
        file.addVariable('dt', dt)
        file.addVariable('inclination', inclination)
        file.addVariable('W', W)
        file.addVariable('S', S)
        file.addVariable('PReqMin', PReqMin)
        file.addVariable('PReqMax', PReqMax)
        file.addVariable('severity', severity)
        file.addVariable('groundCovered', groundCovered)
        file.addVariable('solarGroundCovered', solarGroundCovered)
        
        finalFilename = 'stitched_' + str(round(preDiveHeight, 1)) + "at" + \
            str(round(postAscentVHor, 1)) + 'to' + str(round(postDiveHeight, 1)) + \
            'at' + str(round(postDiveVHor, 1)) + '.dat'
                  
        file.save(finalFilename)

    if plotResult:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(solarGroundCovered / 1e3, heightTotal / 1e3, 'g', label='solar distance', linewidth=2)
        ax.plot(groundCovered / 1e3, heightTotal / 1e3, 'b', label='venusian distance')
        ax.set_xlabel('Ground track [km]')
        ax.set_ylabel('Height [km]')
        ax.legend()
        ax.grid(True)
    
        prevTime = int(timeTotal[0])
    
        for i in range(1, len(timeTotal)):
            curTime = int(timeTotal[i])
    
            if prevTime % (30 * 60) > 15 * 60 and curTime % (30 * 60) < 15 * 60:
                ax.text(groundCovered[i] / 1e3, heightTotal[i] / 1e3,
                    TrackCommon.FormatTime(curTime, 'hh:mm'))
    
            prevTime = curTime
        
    return finalFilename

def DetermineAreaInitialGuess(time, height, alpha, gamma, power, latitude,
                              longitude, powerEfficiency, PRequired, areaRatio,
                              settings):
    # Determine the lower bound and the upper band of the required area
    # estimations.
    #print(' * minimum height:', np.min(height))
    #print(' * maximum height:', np.max(height))
    effDir, effInd, effUni = powerEfficiency.getPowerEfficiency(np.max(height),
        latitude, longitude, 0, 0)

    areaLowerBound = settings.PConsumption / (effDir + effInd + effUni * areaRatio) / settings.fluxVenus

    effDir, effInd, effUni = powerEfficiency.getPowerEfficiency(np.min(height),
        latitude, longitude, 0, 0)

    areaUpperBound = (settings.PConsumption + PRequired) / \
        ((effDir + effInd + effUni * areaRatio) * settings.efficiencyCharging *
        settings.efficiencyDischarging * settings.efficiencyPower) / settings.fluxVenus

    return areaLowerBound, areaUpperBound

def DetermineArea(time, height, alpha, gamma, power, latitude, longitude,
                  areaRatio, PRequired, powerEfficiency, settings, relax=0.8,
                  plotResults=True):
    print(TrackCommon.StringHeader("Estimating solar area and capacity", 60))
    
    lowerBound, upperBound = DetermineAreaInitialGuess(time, height, alpha,
        gamma, power, latitude, longitude, powerEfficiency, PRequired, areaRatio,
        settings)

    center = (lowerBound + upperBound) / 2.0
    capacity = np.zeros(center.shape)
    numBisections = 64
    numSubUpdate = 2
    numUpdate = 16

    #print(' * lower bounds:', lowerBound)
    #print(' * upper bounds:', upperBound)
    #print(' * initial center:', center)

    # Predetermine efficiencies as they don't change for a given dataset
    effDir = np.zeros([len(time)])
    effInd = np.zeros([len(time)])
    effUni = np.zeros([len(time)])

    for i in range(0, len(time)):
        effDir[i], effInd[i], effUni[i] = powerEfficiency.getPowerEfficiency(
            height[i], latitude, longitude, alpha[i] / 180.0 * np.pi, gamma[i])

    PTotalRequired = power / settings.efficiencyPropellers + settings.PConsumption

    # Predetermine the efficiency multiplication factor for each of the area
    # ratios
    PAvailableFactor = np.zeros([len(areaRatio), len(time)])

    for i in range(0, len(areaRatio)):
        PAvailableFactor[i] = (effDir + effInd + effUni * areaRatio[i])

    estimator = TimeEstimator.TimeEstimator(numBisections)
    estimator.startTiming()

    for iBisection in range(0, numBisections):
        # Loop through all provided area ratios
        #print('Iteration:', iBisection)
        estimator.startIteration(iBisection)

        for iRatio in range(0, len(areaRatio)):
            # For the current area ratio determine the points where the power
            # available changes sign with respect to the power required
            PTotalRatio = (effDir + effInd + effUni * areaRatio[iRatio]) * \
                settings.efficiencyPower * settings.fluxVenus
            PTotalAvailable = center[iRatio] * PTotalRatio

            #print(' > Area ratio', iRatio, ':', areaRatio[iRatio])

            # Decompose the array of boolean values into ranges
            rangeExcess = []
            rangeLack = []

            currentExcess = PTotalAvailable[0] > PTotalRequired[0]
            currentIndex = 0

            for iTime in range(1, len(time)):
                if (PTotalAvailable[iTime] > PTotalRequired[iTime]) != currentExcess:
                    if currentExcess:
                        rangeExcess.append([currentIndex, iTime])
                    else:
                        rangeLack.append([currentIndex, iTime])

                    currentExcess = (PTotalAvailable[iTime] > PTotalRequired[iTime])
                    currentIndex = iTime

            if currentExcess:
                rangeExcess.append([currentIndex, len(time)])
            else:
                rangeLack.append([currentIndex, len(time)])

            # Use the ranges to perform the required integrations
            GeneratedExcess = 0.0
            GeneratedLack = 0.0
            RequiredExcess = 0.0
            RequiredLack = 0.0

            #print(' > * Found excess:', rangeExcess)
            #print(' > * Found lack:  ', rangeLack)

            for curExcess in rangeExcess:
                GeneratedExcess += scp_int.simps(
                    PTotalRatio[curExcess[0]:curExcess[1]],
                    time[curExcess[0]:curExcess[1]])
                RequiredExcess += scp_int.simps(
                    PTotalRequired[curExcess[0]:curExcess[1]],
                    time[curExcess[0]:curExcess[1]])

            for curLack in rangeLack:
                GeneratedLack += scp_int.simps(
                    PTotalRatio[curLack[0]:curLack[1]],
                    time[curLack[0]:curLack[1]])
                RequiredLack += scp_int.simps(
                    PTotalRequired[curLack[0]:curLack[1]],
                    time[curLack[0]:curLack[1]])

            #print(' > * Generated excess:', GeneratedExcess)
            #print(' > * Required excess: ', RequiredExcess)
            #print(' > * Generated lack:', GeneratedLack)
            #print(' > * Required lack:', RequiredLack)

            # Use the integral values to determine the area
            Atop = (RequiredLack + settings.efficiencyCharging * settings.efficiencyDischarging *
                RequiredExcess) / (settings.efficiencyPower * (settings.efficiencyDischarging *
                settings.efficiencyCharging * GeneratedExcess + GeneratedLack))

            # Choose new area based on calculated value
            if Atop > center[iRatio]:
                # Increase area
                lowerBound[iRatio] = center[iRatio]
            else:
                # Decrease area
                upperBound[iRatio] = center[iRatio]

            #print(' > * Guess:', round(center[iRatio], 5),
            #    ', gave:', round(Atop, 5), ", new:",
            #    round((lowerBound[iRatio] + upperBound[iRatio]) / 2.0, 5))

            # Calculate the associated battery capacity
            capacity[iRatio] = settings.efficiencyCharging * (settings.efficiencyPower *
                GeneratedExcess * center[iRatio] - RequiredExcess)
            #print(' > * Capacity:', round(capacity[iRatio] / 1e6, 4), 'MJ')

        # Calculate new centers
        centerOld = center
        center = (lowerBound + upperBound) / 2.0
        center = centerOld + (center - centerOld) * relax

        estimator.finishedIteration(iBisection)
        
        if (iBisection + 1) % numSubUpdate == 0:
            print('.', end='')
            
        if (iBisection + 1) % numUpdate == 0:
            print(' spent:', estimator.getTotalElapsed(),
                  ', remaining:', estimator.getEstimatedRemaining())


    return center, capacity

def GenerateThrustFile(filename, thrustFilename):
    atm = Atmosphere.Atmosphere()
    file = TrackStorage.DataStorage()
    file.load(filename)

    # Generate dynamic pressure and density values
    time = file.getVariable('timeTotal').getValues()
    height = file.getVariable('heightTotal').getValues()
    vInf = file.getVariable('vInfTotal').getValues()
    latitude = file.getVariable('latitude').getValues()
    longitude = file.getVariable('longitude').getValues()
    severity = file.getVariable('severity').getValues()
    power = file.getVariable('powerTotal').getValues()
    speedOfSound = atm.speedOfSound(height, latitude, longitude)
    temp = TrackCommon.AdjustSeverity(atm.temperature(height, latitude, longitude), severity)
    density = TrackCommon.AdjustSeverity(atm.density(height, latitude, longitude), severity)
    qInf = 0.5 * density * np.power(vInf, 2.0)

    np.savetxt(thrustFilename, np.transpose([time, height, power, vInf,
        density, qInf, speedOfSound, temp]), '%10.8f', delimiter=';',
        header='time [s]; height [m]; power [W]; vInf [m/s]; density [kg/m3]; ' +
        'qInf [Pa]; a [m/s]; T [K]')

def AnalyzePower(time, height, vHor, vVer, vInf, alpha, gamma, power, latitude,
                 longitude, settings, atmosphere):
    # Go through all segments and seek out the points where there is too little
    # power available to power all subsystems directly
    powerEfficiency = TrackPower.TrackPower(settings, atmosphere)
    powerAvailable = np.zeros([len(time)])

    for i in range(0, len(time)):
        powerAvailable = powerEfficiency.getPowerEfficiency(height[i], latitude,
            longitude, alpha[i], gamma[i])

def TestDetermineArea():
    settings = TrackSettings.Settings()
    atmosphere = Atmosphere.Atmosphere()
    powerEfficiency = TrackPower.TrackPower(settings, atmosphere)
    file = TrackStorage.DataStorage()
    file.load('stitched_62000.0at7.8to38000.0at-25.0.dat')
    time = file.getVariable('timeTotal').getValues()
    height = file.getVariable('heightTotal').getValues()
    vHor = file.getVariable('vHorTotal').getValues()
    vVer = file.getVariable('vVerTotal').getValues()
    alpha = file.getVariable('alphaTotal').getValues()
    gamma = file.getVariable('gammaTotal').getValues()
    power = file.getVariable('powerTotal').getValues()
    latitude = file.getVariable('latitude').getValues()
    longitude = file.getVariable('longitude').getValues()
    PRequired = file.getVariable('PReqMax').getValues()
    areaRatio = np.linspace(2, 2.8, 9)

    areas, capacity = DetermineArea(time, height, alpha, gamma, power, latitude,
        longitude, areaRatio, PRequired, powerEfficiency, settings)

    for i in range(0, len(areaRatio)):
        print('Area ratio:', areaRatio[i])
        print(' > Atop:   ', round(areas[i], 2), 'm2')
        print(' > Abottom:', round(areas[i], 2), 'm2')
        print(' > Atotal: ', round(areaRatio[i] * areas[i], 2), 'm2')
        print(' > capacity:', round(capacity[i] / 1e6, 3), 'MJ')

        batWeight, batVolume = AircraftBatteries.AircraftBatterySizing(
            capacity[i] / 3600, settings.batteryDepthOfDischarge,
            settings.batterySafetyFactor)

        solarPanelWeight = areaRatio[i] * areas[i] * settings.specificWeightPanels
        print(' > battery weight:', round(batWeight, 3), 'kg')
        print(' > solar panel weight:', round(solarPanelWeight, 3), 'kg')
        print(' > total weight:', round(solarPanelWeight + batWeight, 3), 'kg')
#StitchTracks(62e3, 3.5, 38e3, -25.0, 10, -5.0, 7.8,
#             0e3, 32e3, 0.20, TrackSettings.Settings(), 0.0)
#TestDetermineArea()
#GenerateThrustFile('stitched_62000.0at7.8to38000.0at-25.0.dat',
#                   'thrust_62000.0at7.8to38000at-25.0.csv')
