# -*- coding: utf-8 -*-
"""
Created on Mon May  9 19:52:15 2016

@author: MaxHenger
"""

import utility as util
import numpy as np
import scipy.integrate as scp_int
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Move to main folder for extermal imports
import os
curDir = os.getcwd()
os.chdir("../../")

import Atmosphere

os.chdir(curDir)

def calculateAveragePowerOverArea(settings, atm, alpha1, h, t, latitude, Veff, etaTot):
    dadt = ( settings.venusOmega * (settings.venusRadius + h) * np.cos(latitude) + \
        Veff) / (settings.venusRadius + h)
    
    if (alpha1 is None):
        alpha1 = -dadt * t / 2.0
        alpha2 = -alpha1
    else:
        alpha2 = alpha1 + dadt * t
        
    etaAtm = atm.solarEfficiency(h, 0.0, 0.0, False)
    
    '''print('h =', h,
          ', alpha1 =', round(alpha1 * 180.0 / np.pi, 3), 
          ', alpha2 =', round(alpha2 * 180.0 / np.pi, 3),
          ', lat =', round(latitude * 180.0 / np.pi),
          ', dadt =', dadt)'''
    
    return settings.venusFlux * etaAtm * etaTot / (t * np.cos(latitude)) * \
        ((settings.venusRadius + h) / \
        (settings.venusOmega * (settings.venusRadius + h) * np.cos(latitude) + Veff) * \
        (np.sin(alpha2) - np.sin(alpha1)) - t * np.power(np.sin(latitude), 2.0))
        
def calculateMissionBounds(settings, atm, 
                           heightUpper, tUpper, POverSUpper, VUpper,
                           heightLower, tLower, POverSLower, VLower,
                           latitudes, eta):
    longitudes = np.zeros(latitudes.shape)
    
    # Precalculate the average required POverS over the upper and lower trakcs
    POverSReq = (tUpper * POverSUpper + tLower * POverSLower) / (tUpper + tLower)
    
    for iLat in range(0, len(latitudes)):
        low = -60.0 * np.pi / 180.0
        high = 60.0 * np.pi / 180.0
            
        for iIter in range(0, 15):
            center = (low + high) / 2.0
            #print('lat =', latitudes[iLat] / np.pi  * 180.0, ', center =', center / np.pi * 180.0)
            
            curUpper = calculateAveragePowerOverArea(settings, atm, center, heightUpper,
                                                        tUpper, latitudes[iLat],
                                                        VUpper, eta)
            curLower = calculateAveragePowerOverArea(settings, atm, center, heightLower,
                                                        tLower, latitudes[iLat],
                                                        VLower, eta)
            
            #print('res up =', curUpper, 'res low =', curLower)
            POverSCur = (curUpper * tUpper + curLower * tLower) / (tUpper + tLower)
            
            if POverSCur > POverSReq:
                '''print('* POverS =', POverSCur, '>', POverSReq,
                      ', lat =', latitudes[iLat] / np.pi * 180.0,
                      ', center =', center / np.pi * 180.0)'''
                low = center
            else:
                '''print('* POverS =', POverSCur, '<', POverSReq,
                      ', lat =', latitudes[iLat] / np.pi * 180.0,
                      ', center =', center / np.pi * 180.0)'''
                high = center
                
        center = (low + high / 2)
        
        if center > 0.0:
            longitudes[iLat] = center
        else:
            longitudes[iLat] = 0.0
        
    return longitudes
                
        
'''def calculateDeltaVLower(settings, atm, longitude, latitude):
    # Retrieve worst case wind speeds (very low at lower lap, very high at
    # upper lap)
    VWindUpper = atm.velocityZonal(settings.heightUpper, latitude, longitude)
    VWindLower = atm.velocityZonal(settings.heightLower, latitude, longitude)
    VWindUpper = np.max(VWindUpper)
    VWindLower = np.min(VWindLower)
    
    # Calculate required upper deltaV
    return settings.ratioTime * (settings.venusRadius + settings.heightLower) / \
        (settings.venusRadius + settings.heightUpper) * \
        (settings.venusOmega * (settings.venusRadius + settings.heightUpper) * np.cos(latitude) + \
        VWindUpper - settings.deltaVUpper) - \
        settings.venusOmega * (settings.venusRadius + settings.heightLower) * np.cos(latitude) - \
        VWindLower'''

def calculateAerodynamicProperties(settings, atm, m, h, V, lat):
    # Settings
    numIntegrationPoints = 150
    
    # Perform funky integration (elliptical lift distribution with a fuselage
    # degrading performance)
    xInt = np.linspace(-3.0 * settings.radiusHeatshield,
                       3.0 * settings.radiusHeatshield,
                       numIntegrationPoints)
    yInt = np.zeros(xInt.shape)
    
    for i in range(0, len(xInt)):
        factor = 1.0
        
        if abs(xInt[i]) < settings.fuselageRadius:
            localX = abs(xInt[i])                
            factorTerm = localX / settings.fuselageRadius
            factor = 1 - (2.0 * factorTerm**3.0 - 3.0 * factorTerm**2.0 + 1.0)**2.0 * 0.5
            
        yInt[i] = factor * np.sqrt(1 - (xInt[i] / (3.0 * settings.radiusHeatshield))**2.0)
        
    integration = scp_int.simps(yInt, xInt)
    
    density = atm.density(h, lat, 0)[0]
    #print('h =', round(h / 1e3, 3), 'km, density =', density, 'kg/m3')
    q = 0.5 * density * V**2.0
    chord = m * settings.venusMu / (settings.venusRadius + h)**2.0 / \
        (q * settings.designLiftCoefficient * integration)
    area = 6.0 * settings.radiusHeatshield * chord
    liftCoefficient = (m * settings.venusMu / (settings.venusRadius + h)**2.0) / \
        (q * area)
    dragCoefficient = settings.minimumDrag + liftCoefficient**2.0 / (np.pi * 
        settings.oswaldFactor * area / chord**2.0)
    
    return [chord, area, liftCoefficient, liftCoefficient * q * area,
            dragCoefficient, dragCoefficient * q * area]
    
def analyzeMissionRange(settings, printResults=False):
    # Some script settings and often used variables
    numVelocityPoints = 50
    atm = Atmosphere.Atmosphere()
    
    # Perform an initial total power and mass estimation
    totalMass = settings.massPayload * settings.ratioMass0 + settings.ratioMass1
    totalPower = settings.powerPayload * settings.ratioPower0 + settings.ratioPower1
    
    # Determine atmospherically related variables
    tLower = settings.timeLower
    tUpper = settings.timeLower * settings.ratioTime
    
    latTicks = np.linspace(0, settings.deltaLatitude, numVelocityPoints)
    lonTicks = np.linspace(-settings.deltaLatitude, settings.deltaLatitude, numVelocityPoints)
    
    lat = np.zeros([numVelocityPoints * numVelocityPoints])
    lon = np.zeros(lat.shape)
    
    for iLat in range(0, len(latTicks)):
        for iLon in range(0, len(lonTicks)):
            i = iLat * numVelocityPoints + iLon
            lat[i] = latTicks[iLat]
            lon[i] = lonTicks[iLon]
    
    VWindLower = atm.velocityZonal(settings.heightLower, lat, lon)
    VWindLower = np.max(VWindLower)
    VWindUpper = atm.velocityZonal(settings.heightUpper, lat, lon)
    VWindUpper = np.max(VWindUpper)
    
    VUpper = VWindUpper - settings.deltaVUpper
    VLower = - (settings.venusRadius + settings.heightLower) * ( \
            (tUpper / tLower) * ( \
                (VUpper - VWindUpper) / (settings.venusRadius + settings.heightUpper) -
                settings.venusOmega * np.cos(settings.deltaLatitude)
            ) - settings.venusOmega * np.cos(settings.deltaLatitude)
        ) + VWindLower
    deltaVLower = VWindLower - VLower
        
    # Determine the flyable range over Venus with the current settings
    deltaAlpha = (settings.venusOmega * (settings.venusRadius + settings.heightUpper) * \
        np.cos(settings.deltaLatitude) + VUpper) / (2 * (settings.venusRadius + settings.heightUpper))

    deltaLongitude = np.arccos(np.cos(deltaAlpha) - np.power(np.sin(settings.deltaLatitude), 2.0) / \
        np.power(np.cos(settings.deltaLatitude), 2.0))
        
    # Use the low er and upper atmospheric solar cell efficiency to make an 
    # initial estimate of surface area and battery capacity
    POverSUpper = calculateAveragePowerOverArea(settings, atm, None, 
                                                settings.heightUpper, tUpper, 
                                                settings.deltaLatitude, VUpper - VWindUpper,
                                                settings.efficiencySolarPanel)
    POverSLower = calculateAveragePowerOverArea(settings, atm, None,
                                                settings.heightLower, tLower,
                                                settings.deltaLatitude, VLower - VWindLower,
                                                settings.efficiencySolarPanel)
    
    boundsLat = np.linspace(0, 75.0 * np.pi / 180.0, 50)
    bounds = calculateMissionBounds(settings, atm, settings.heightUpper,
                                    tUpper, POverSUpper, VUpper, settings.heightLower,
                                    tLower, POverSLower, VLower, boundsLat,
                                    settings.efficiencySolarPanel)
    
    bounds = np.asarray([boundsLat, bounds])
    
    # Perform aerodynamic estimation
    aeroUpper = calculateAerodynamicProperties(settings, atm, totalMass, settings.heightUpper, 
                                               VUpper, settings.deltaLatitude)
    aeroLower = calculateAerodynamicProperties(settings, atm, totalMass, settings.heightLower,
                                               VLower, settings.deltaLatitude)    
    
    propPowerUpper = 0.5 * aeroUpper[5] * abs(VUpper) * \
        (np.sqrt(aeroUpper[1] * aeroUpper[4] / settings.propellorArea + 1.0) + 1.0)
    propPowerLower = 0.5 * aeroLower[5] * abs(VLower) * \
        (np.sqrt(aeroLower[1] * aeroLower[4] / settings.propellorArea + 1.0) + 1.0)
    
    SSolarCell = (((totalPower + propPowerUpper) * tUpper) / settings.efficiencyPower + 
        ((totalPower + propPowerLower) * tLower) / (settings.efficiencyPower * settings.efficiencyCharging)) / \
        (POverSUpper * tUpper + POverSLower * tLower / settings.efficiencyCharging)
        
    CBattery = (totalPower + propPowerLower) * tLower / settings.efficiencyPower - \
        POverSLower * SSolarCell * tLower
        
    VGUpper = (settings.venusRadius / (settings.venusRadius + settings.heightUpper)) * \
        (VUpper - VWindUpper) - settings.venusRadius * settings.venusOmega * np.cos(settings.deltaLatitude)
    VGLower = (settings.venusRadius / (settings.venusRadius + settings.heightLower)) * \
        (VLower - VWindLower) - settings.venusRadius * settings.venusOmega * np.cos(settings.deltaLatitude)
    
    if printResults == True:
        print('----- general parameters -----')
        print('total mass estimate =', totalMass, 'kg (payload =', settings.massPayload, 'kg)')
        print('total power estimate =', totalPower, 'W (payload =', settings.powerPayload, 'W)')
        
        print('----- temporal parameters -----')
        print('time upper =', tUpper, 's')
        print('time lower =', tLower, 's')
        print('velocity upper =', round(VWindUpper, 3), 
              '-', round(settings.deltaVUpper, 3), 
              '=', round(VUpper, 3), 'm/s')
        print('velocity lower =', round(VWindLower, 3), 
              '-', round(deltaVLower, 3), 
              '=', round(VLower, 3), 'm/s')
        print('ground velocity upper =', round(VGUpper, 3), 'm/s')
        print('ground velocity lower =', round(VGLower, 3), 'm/s')
        
        print('----- aerodynamic parameters -----')
        print('chord lower =', round(aeroLower[0], 3), 'm')
        print('area lower =', round(aeroLower[1], 3), 'm2')
        print('CL lower =', round(aeroLower[2], 3))
        print('Lift lower =', round(aeroLower[3] / 1e3, 3), 'kN')
        print('CD lower =', round(aeroLower[4], 3))
        print('Drag lower =', round(aeroLower[5] / 1e3, 3), 'kN')
        print('chord upper =', round(aeroUpper[0], 3), 'm')
        print('area upper =', round(aeroUpper[1], 3), 'm2')
        print('CL upper =', round(aeroUpper[2], 3))
        print('Lift upper =', round(aeroUpper[3] /1e3, 3), 'kN')
        print('CD upper =', round(aeroUpper[4], 3))
        print('Drag upper =', round(aeroUpper[5] /1e3, 3), 'kN')
        
        print('----- power parameters -----')
        print('upper P/S =', POverSUpper, 'W/m2')
        print('lower P/S =', POverSLower, 'W/m2')
        print('solar cell area =', SSolarCell, 'm2')
        print('upper prop power =', round(propPowerUpper /1e3, 3), 'kW')
        print('lower prop power =', round(propPowerLower / 1e3, 3), 'kW')
        
        print('battery =', CBattery / 1e6, 'MJ =', CBattery / settings.specificCapacityBattery, 'kg')
        print('latitude range =', settings.deltaLatitude * 180.0 / np.pi, 'deg')
        print('longitude range (at highest latitude) =', deltaLongitude * 180.0 / np.pi, 'deg')
        
    return {
        'totalMass': totalMass,
        'totalPower': totalPower,
        'tUpper': tUpper,
        'VWindUpper': VWindUpper,
        'deltaVUpper': settings.deltaVUpper,
        'VUpper': VUpper,
        'VGUpper': VGUpper,
        'POverSUpper': POverSUpper,
        'chordUpper': aeroUpper[0],
        'areaUpper': aeroUpper[1],
        'tLower': tLower,
        'VWindLower': VWindLower,
        'deltaVLower': deltaVLower,
        'VLower': VLower,
        'VGLower': VGLower,
        'POverSLower': POverSLower,
        'chordLower': aeroLower[0],
        'areaLower': aeroLower[1],
        'SSolarCell': SSolarCell,
        'CBattery': CBattery,
        'bounds': bounds,
        'heightUpper': settings.heightUpper,
        'heightLower': settings.heightLower
    }
            
def graphMissionRange():
    settings = util.settings()
    
    deltaLatitude = np.linspace(5.0, 45.0, 5) / 180.0 * np.pi
    results = []
    
    for iLat in range(0, len(deltaLatitude)):
        settings.deltaLatitude = deltaLatitude[iLat]
        results.append(analyzeMissionRange(settings))
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    util.drawVenus(ax, np.linspace(-90, 90, 19), np.linspace(-180, 180, 19),
                   './data/venus.png')
    
    lines = []
    text = []
    for iLat in range(0, len(deltaLatitude)):
        iNonZero = 0
        
        while iNonZero < len(results[iLat][-1][1]):
            if abs(results[iLat][-1][1][iNonZero]) < 1e-8:
                break
            
            iNonZero += 1
        
        resLat = np.zeros(iNonZero * 4 + 1)
        resLon = np.zeros(iNonZero * 4 + 1)
        
        print('lat =', deltaLatitude[iLat] * 180.0 / np.pi,
              ', nonzero =', iNonZero, ', len =', len(resLat))
        
        for i in range(0, iNonZero):
            resLat[0 * iNonZero + i] = results[iLat][-1][0][i]
            resLat[1 * iNonZero + i] = results[iLat][-1][0][iNonZero - 1 - i]
            resLat[2 * iNonZero + i] = -results[iLat][-1][0][i]
            resLat[3 * iNonZero + i] = -results[iLat][-1][0][iNonZero - 1 - i]
            resLon[0 * iNonZero + i] = results[iLat][-1][1][i]
            resLon[1 * iNonZero + i] = -results[iLat][-1][1][iNonZero - 1 - i]
            resLon[2 * iNonZero + i] = -results[iLat][-1][1][i]
            resLon[3 * iNonZero + i] = results[iLat][-1][1][iNonZero - 1 - i]
            
        resLat[-1] = resLat[0]
        resLon[-1] = resLon[0]
        
        line, = ax.plot(resLat * 180.0 / np.pi, resLon * 180.0 / np.pi)
        lines.append(line)
        text.append('lat =' + str(round(deltaLatitude[iLat] * 180.0 / np.pi)))
            
    ax.legend(lines, text)
    
def analyzeSingleMission():
    analyzeMissionRange(util.settings(), True)
    
def analyzeMissionParameters(filename):
    numPoints = 30
    higherAltitude = np.linspace(30000, 70000, numPoints)
    lowerAltitude = np.linspace(30000, 70000, numPoints)
    #timeRatio = np.linspace(2.0, 8.0, 7)
    timeRatio = [1.5, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
    deltaVUpper = np.linspace(-50, 50, 11)
    
    settings = util.settings()
    
    iLowestAltLow = 0
    iLowestAltHigh = 0
    iLowestTime = 0
    iLowestV = 0
    lowestS = 1e9
    
    allResults = []
    
    for iAltLow in range(0, len(lowerAltitude)):
        for iAltHigh in range(0, len(higherAltitude)):
            print('processing', iAltLow * len(higherAltitude) + iAltHigh + 1, 
                  'of', len(lowerAltitude) * len(higherAltitude))
            
            for iTime in range(0, len(timeRatio)):
                for iV in range(0, len(deltaVUpper)):
                    settings.heightLower = lowerAltitude[iAltLow]
                    settings.heightUpper = higherAltitude[iAltHigh]
                    settings.ratioTime = timeRatio[iTime]
                    settings.deltaVUpper = deltaVUpper[iV]
                    results = analyzeMissionRange(settings)
                    
                    allResults.append(results)
                    
                    #if results[-3] < 0:
                        #analyzeMissionRange(settings, True)
                        #raise ValueError('INVALID VALUES')
                        
                    if results['SSolarCell'] > 0 and results['SSolarCell'] < lowestS and \
                        results['chordLower'] < 3.5 and results['chordUpper'] < 3.5 and \
                        results['VUpper'] < 100.0 and results['VLower'] < 100.0:
                        print('new area =', round(results['SSolarCell'], 3), 'm2',
                              ', alt_u =', round(higherAltitude[iAltHigh] / 1e3, 3), 'km',
                              ', alt_l =', round(lowerAltitude[iAltLow] / 1e3, 3), 'km',
                              ', f_time =', round(timeRatio[iTime] * 1e2, 3), '%',
                              ', deltaV =', round(deltaVUpper[iV], 3), 'm/s')
                        lowestS = results['SSolarCell']
                        iLowestAltLow = iAltLow
                        iLowestAltHigh = iAltHigh
                        iLowestTime = iTime
                        iLowestV = iV
    
    settings.heightLower = lowerAltitude[iLowestAltLow]
    settings.heightUpper = higherAltitude[iLowestAltHigh]
    settings.ratioTime = timeRatio[iLowestTime]
    settings.deltaVUpper = deltaVUpper[iLowestV]
    analyzeMissionRange(settings, True)
    
    dictionary = util.dictionaryIO()
    dictionary.setAxes(['heightLower', 'heightUpper', 'ratioTime', 'deltaVUpper'],
                       [lowerAltitude, higherAltitude, timeRatio, deltaVUpper])
    dictionary.setResults(allResults)
    dictionary.save(filename)
    
# plotResults is the function that loads the information saved with the
# analyzeMissionParameters function and plots the results graphically.
# Input:
#   ax1: the index of the axis to plot along the x-axis
#   ax2: the index of the axis to plot along the y-axis
#   name: name of the variable to plot
#   filename: the name of the file in which all results are saved
#   minimum: an array of 3-tuples: (variable name, contour to plot, 
#       line color, line specification)
#   maximum: an array of 3-tuples
#   lowest: when true the variable indicated by 'name' will be minimized. If
#       false then the variable indicated by 'name' will be maximized
def plotResults(ax1, ax2, name, filename, minimum=[], maximum=[], lowest=True):
    dictionary = util.dictionaryIO()
    dictionary.load(filename)
    
    axis1Name = dictionary.getAxisName(ax1)
    axis1 = dictionary.getAxisValues(ax1)
    axis2Name = dictionary.getAxisName(ax2)
    axis2 = dictionary.getAxisValues(ax2)
    
    remaining = [0, 1, 2, 3]
    if ax1 < ax2:
        del(remaining[ax2])
        del(remaining[ax1])
    else:
        del(remaining[ax1])
        del(remaining[ax2])
    
    remainingAxis1 = dictionary.getAxisValues(remaining[0])
    remainingAxis1Name = dictionary.getAxisName(remaining[0])
    remainingAxis2 = dictionary.getAxisValues(remaining[1])
    remainingAxis2Name = dictionary.getAxisName(remaining[1])
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    results = np.zeros([len(axis2), len(axis1)])
    minResults = np.zeros([len(minimum), len(axis2), len(axis1)])
    maxResults = np.zeros([len(maximum), len(axis2), len(axis1)])
    
    iAbsBestAxis1 = 0
    iAbsBestAxis2 = 0
    iAbsBestRem1 = 0
    iAbsBestRem2 = 0
    absBestMetric = 1e9
    
    if lowest == False:
        metric = -metric
    
    for iAxis1 in range(0, len(axis1)):
        for iAxis2 in range(0, len(axis2)):
            # Search for the highest/lowest metric in the remaining variables
            metric = 1e9
            iBestRem1 = 0
            iBestRem2 = 0
            
            if lowest == False:
                metric = -metric
                
            for iRem1 in range(0, len(remainingAxis1)):
                for iRem2 in range(0, len(remainingAxis2)):
                    indices = [0, 0, 0, 0]
                    indices[ax1] = iAxis1
                    indices[ax2] = iAxis2
                    indices[remaining[0]] = iRem1
                    indices[remaining[1]] = iRem2
                    
                    local = dictionary.getValue(indices)
                    
                    if lowest == True:
                        if local[name] < metric:
                            metric = local[name]
                            iBestRem1 = iRem1
                            iBestRem2 = iRem2
                            
                            # Check if it is the absolute best metric
                            if metric < absBestMetric:
                                valid = True
                                for iMin in range(0, len(minimum)):
                                    if local[minimum[iMin][0]] < minimum[iMin][1]:
                                        valid = False
                                        break
                                    
                                if valid == True:
                                    for iMax in range(0, len(maximum)):
                                        if local[maximum[iMax][0]] > maximum[iMax][1]:
                                            valid = False
                                            break
                                        
                                if valid == True:
                                    iAbsBestAxis1 = iAxis1
                                    iAbsBestAxis2 = iAxis2
                                    iAbsBestRem1 = iRem1
                                    iAbsBestRem2 = iRem2
                                    absBestMetric = metric
                                    
                            
                    else:
                        if local[name] > metric:
                            metric = local[name]
                            iBestRem1 = iRem1
                            iBestRem2 = iRem2
                            
            indices = [0, 0, 0, 0]
            indices[ax1] = iAxis1
            indices[ax2] = iAxis2
            indices[remaining[0]] = iBestRem1
            indices[remaining[1]] = iBestRem2
            local = dictionary.getValue(indices)
            
            for iMin in range(0, len(minimum)):
                minResults[iMin, iAxis2, iAxis1] = local[minimum[iMin][0]]
            
            for iMax in range(0, len(maximum)):
                maxResults[iMax, iAxis2, iAxis1] = local[maximum[iMax][0]]
                            
            results[iAxis2, iAxis1] = metric
            
    img = ax.imshow(results, extent=[min(axis1), max(axis1), max(axis2), min(axis2)], aspect='auto')
    fig.colorbar(img)
    
    #print(maxResults[0])
    
    # Plot minimum and maximum contours
    for i in range(0, len(minimum)):
        ax.contourf(axis1, axis2, minResults[i], levels=[-1e9, minimum[i][1]],
                    hatches=['//'], colors=None, fill=False)
        ax.contour(axis1, axis2, minResults[i], levels=[minimum[i][1]],
                   colors=minimum[i][2], linestyles=minimum[i][3])
    
    for i in range(0, len(maximum)):
        print
        ax.contourf(axis1, axis2, maxResults[i], levels=[maximum[i][1], 1e9],
                    hatches=['//'], colors=None, fill=False)
        ax.contour(axis1, axis2, maxResults[i], levels=[maximum[i][1]],
                   colors=maximum[i][2], linestyles=maximum[i][3])
        
    # Plot the design point
    ax.plot(axis1[iAbsBestAxis1], axis2[iAbsBestAxis2], 'ro', markersize=30)
        
    # Print results of absolute best
    indices = [0, 0, 0, 0]
    indices[ax1] = iAbsBestAxis1
    indices[ax2] = iAbsBestAxis2
    indices[remaining[0]] = iAbsBestRem1
    indices[remaining[1]] = iAbsBestRem2
    local = dictionary.getValue(indices)
    print(' ****** Results:')
    print(' *** Parameters:')
    print(axis1Name, '=', axis1[iAbsBestAxis1])
    print(axis2Name, '=', axis2[iAbsBestAxis2])
    print(remainingAxis1Name, '=', remainingAxis1[iAbsBestRem1])
    print(remainingAxis2Name, '=', remainingAxis2[iAbsBestRem2])
    print(' *** Temporal values:')
    print('tUpper =', round(local['tUpper'], 3), 's')
    print('tLower =', round(local['tLower'], 3), 's')
    print(' *** Velocity values:')
    print('VUpper =', round(local['VUpper'], 3), 'm/s')
    print('VLower =', round(local['VLower'], 3), 'm/s')
    print('VGUpper =', round(local['VGUpper'], 3), 'm/s')
    print('VGLower =', round(local['VGLower'], 3), 'm/s')
    print(' *** Geometric values:')
    print('heightUpper =', round(local['heightUpper'] / 1e3, 3), 'km')
    print('heightLower =', round(local['heightLower'] / 1e3, 3), 'km')
    print('chordUpper =', round(local['chordUpper'], 3), 'm')
    print('chordLower =', round(local['chordLower'], 3), 'm')
    print('areaUpper =', round(local['areaUpper'], 3), 'm2')
    print('areaLower =', round(local['areaLower'], 3), 'm2')
    print(' *** Power values:')
    print('SSolarCell =', round(local['SSolarCell'], 3), 'm2')
    print('POverSUpper =', round(local['POverSUpper'], 5), 'W/m2')
    print('POverSLower =', round(local['POverSLower'], 5), 'W/m2')
    print('CBattery =', round(local['CBattery'] / 1e6, 5), 'MJ')
    print('mBattery =', round(local['CBattery'] / 0.46e6, 5), 'kg')
    ax.set_xlabel(axis1Name)
    ax.set_ylabel(axis2Name)
    ax.set_title(name)
    
#analyzeMissionRange(util.settings())

def __testLatitudeLongitudeLocal__():
    latitude = np.linspace(0 * np.pi / 180.0, 80.0 * np.pi / 180, 40)
    longitude = np.linspace(10 * np.pi / 180.0, 80.0 * np.pi / 180.0, 20)
    res = np.zeros([len(latitude), len(longitude)])
    
    for iLat in range(0, len(latitude)):
        for iLon in range(0, len(longitude)):
            precos = np.cos(latitude[iLat])**2.0 * np.sin(longitude[iLon]) + \
                np.sin(latitude[iLat])**2.0
            print('lat=',latitude[iLat],'lon=',longitude[iLon],'precos=',precos)
            res[iLat,iLon] = np.arccos(
                np.cos(latitude[iLat])**2.0 * np.sin(longitude[iLon]) + \
                np.sin(latitude[iLat])**2.0
            )
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    gridLat, gridLon = np.meshgrid(latitude * 180.0 / np.pi, longitude * 180.0 / np.pi)
    ax.plot_wireframe(gridLat, gridLon, np.transpose(res) * 180.0 / np.pi)
    ax.set_xlabel('latitude')
    ax.set_ylabel('longitude')
    ax.grid(True)
 
#analyzeMissionRange(util.settings())
#graphMissionRange()
#analyzeSingleMission()
analyzeMissionParameters("results.txt")
plotResults(2, 3, 'SSolarCell', 'results.txt', 
          minimum=[
              ['chordLower', 0.0, 'g', 'solid'],
              ['chordUpper', 0.0, 'g', 'dashed'],
          ],
          maximum=[
              ['chordLower', 4.0, 'r', 'solid'],
              ['chordUpper', 4.0, 'r', 'dashed'],
              ['VLower', 125.0, 'white', 'solid'],
              ['VUpper', 125.0, 'white', 'dashed'],
          ], lowest=True)
#__testLatitudeLongitudeLocal__()