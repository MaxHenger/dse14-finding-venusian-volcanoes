#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 16:57:17 2016

@author: MaxHenger
"""

import numpy as np
import TrackCommon
import TrackContour
import TrackIO
import TrackLookup
import matplotlib.pyplot as plt

# DetermineLengths is a design tool to figure out the relation between the
# aeroshell design and the wing planform design. For a given required half-wing
# surface area, an aeroshell bottom radius, a maximum aeroshell height and a
# defined root chord the relationship between the mid-chord length and the
# aeroshell height, and the wing planform design to achieve the desired
# required area is displayed. If a solution space exists then the contours will
# indicate it.
#
# Input:
#   - SReq: Half of the required wing planform area in m2 (i.e. the area that
#       is formed by a single half of a set of wings)
#   - rAeroBottom: The aeroshell lower radius in m
#   - hAeroMax: The maximum aeroshell height in m
#   - cRoot: The root chord in m
#   - [desCi:] The design point mid-chord length to print the results for in m
#   - [desHae:] The design point aeroshell height to print the results for in m
def DetermineLengths(SReq, rAeroBottom, hAeroMax, cRoot, desCi=None, desHae=None):
    reductionFactor = 0.9
    numPoints = 200
    thickness = 0.1

    axisCi = np.linspace((1.0 - reductionFactor) * cRoot, cRoot, numPoints)
    axisHae = np.linspace((1 - reductionFactor) * hAeroMax, hAeroMax, numPoints)
    resultsB1 = np.zeros([numPoints, numPoints])
    resultsB2 = np.zeros([numPoints, numPoints])
    resultsSpan = np.zeros([numPoints, numPoints])
    resultsTip = np.zeros([numPoints, numPoints])
    resultsRadiusTop = np.zeros([numPoints, numPoints])
    resultsTaperRatio = np.zeros([numPoints, numPoints])

    for iH in range(0, len(axisHae)):
        for iC in range(0, len(axisCi)):
            D = np.sqrt(np.power((2 * thickness * axisCi[iC]), 2.0) -
                4 * (np.power(thickness * axisCi[iC], 2.0) +
                     0.25 * np.power(axisCi[iC], 2.0) -
                     np.power(rAeroBottom, 2.0)))
            b1 = - thickness * axisCi[iC] + D / 2
            b2 = np.sqrt(np.power(b1 + 0.5 * thickness * axisCi[iC], 2.0) +
                         np.power(axisHae[iH], 2.0))
            Ct = 2 * SReq / b2 - b1 / b2 * (cRoot + axisCi[iC]) - axisCi[iC]
            rae_t = Ct * (thickness + 1.0) / 2.0

            resultsB1[iH, iC] = b1
            resultsB2[iH, iC] = b2
            resultsSpan[iH, iC] = 2 * (b1 + b2)
            resultsTip[iH, iC] = Ct
            resultsRadiusTop[iH, iC] = rae_t
            resultsTaperRatio[iH, iC] = Ct / axisCi[iC]

    fig = plt.figure()
    boundsSpan1L, boundsSpan1R = TrackCommon.ImageAxes(0.05, 0.47, 0.53, 0.97)
    boundsSpan2L, boundsSpan2R = TrackCommon.ImageAxes(0.53, 0.95, 0.53, 0.97)
    boundsTipL, boundsTipR = TrackCommon.ImageAxes(0.05, 0.47, 0.03, 0.47)
    radiusTopL, radiusTopR = TrackCommon.ImageAxes(0.53, 0.95, 0.03, 0.47)

    axSpan1Data = fig.add_axes(boundsSpan1L)
    axSpan1Legend = fig.add_axes(boundsSpan1R)
    TrackCommon.PlotImage(fig, axSpan1Data, axSpan1Legend, axisCi, r'$C_i \; [m]$',
        axisHae, r'$h_{ae} \; [m]$', resultsB1, r'$b_1 \; [m]$')

    axSpan2Data = fig.add_axes(boundsSpan2L)
    axSpan2Legend = fig.add_axes(boundsSpan2R)
    TrackCommon.PlotImage(fig, axSpan2Data, axSpan2Legend, axisCi, r'$C_i \; [m]$',
        axisHae, r'$h_{ae} \; [m]$', resultsB2, r'$b_2 \; [m]$')


    axTipData = fig.add_axes(boundsTipL)
    axTipLegend = fig.add_axes(boundsTipR)
    TrackCommon.PlotImage(fig, axTipData, axTipLegend, axisCi, r'$C_i \; [m]$',
        axisHae, r'$h_{ae} \; [m]$', resultsTaperRatio, r'$c_t \; [m]$')

    axRadiusTopData = fig.add_axes(radiusTopL)
    axRadiusTopLegend = fig.add_axes(radiusTopR)
    TrackCommon.PlotImage(fig, axRadiusTopData, axRadiusTopLegend, axisCi,
        r'$C_i \; [m]$', axisHae, r'$h_{ae} \; [m]$', resultsRadiusTop,
        r'$r_{ae,t} \; [m]$')

    contour = TrackContour.Contour()
    contour.combineData(axisCi, axisHae, [resultsB1, resultsTip,
        resultsRadiusTop], [1.4, 0, 1.0], [3.0, cRoot, 1.1])

    for i in range(0, contour.getNumContours()):
        curContour = contour.getContour(i)
        curIncluder = curContour.getIncluder()

        axSpan1Data.plot(curIncluder[:, 0], curIncluder[:, 1], 'g--')
        axSpan2Data.plot(curIncluder[:, 0], curIncluder[:, 1], 'g--')
        axTipData.plot(curIncluder[:, 0], curIncluder[:, 1], 'g--')
        axRadiusTopData.plot(curIncluder[:, 0], curIncluder[:, 1], 'g--')

    if not (desCi is None) and not (desHae is None):
        iC = TrackCommon.Find1DBisectionAscending(axisCi, desCi)
        iHae = TrackCommon.Find1DBisectionAscending(axisHae, desHae)

        b1 = 0
        b2 = 0
        cTip = 0
        aeroshellTop = 0

        if iC == len(axisCi) - 1:
            if iHae == len(axisHae) - 1:
                # At a corner, cannot interpolate
                b1 = resultsB1[iHae, iC]
                b2 = resultsB2[iHae, iC]
                cTip = resultsTip[iHae, iC]
                aeroshellTop = resultsRadiusTop[iHae, iC]
            else:
                # Interpolate along height axis
                b1 = TrackCommon.Lerp(axisHae[iHae], resultsB1[iHae, iC],
                                      axisHae[iHae + 1], resultsB1[iHae + 1, iC], desHae)
                b2 = TrackCommon.Lerp(axisHae[iHae], resultsB2[iHae, iC],
                                      axisHae[iHae + 1], resultsB2[iHae + 1, iC], desHae)
                cTip = TrackCommon.Lerp(axisHae[iHae], resultsTip[iHae, iC],
                                        axisHae[iHae + 1], resultsTip[iHae + 1, iC], desHae)
                aeroshellTop = TrackCommon.Lerp(axisHae[iHae], resultsRadiusTop[iHae, iC],
                                                axisHae[iHae + 1], resultsRadiusTop[iHae + 1, iC], desHae)
        elif iHae == len(axisHae) - 1:
            # Interpolate along mid-chord axis
            b1 = TrackCommon.Lerp(axisCi[iC], resultsB1[iHae, iC],
                                  axisCi[iC + 1], resultsB1[iHae, iC + 1], desCi)
            b2 = TrackCommon.Lerp(axisCi[iC], resultsB2[iHae, iC],
                                  axisCi[iC + 1], resultsB2[iHae, iC + 1], desCi)
            cTip = TrackCommon.Lerp(axisCi[iC], resultsTip[iHae, iC],
                                    axisCi[iC + 1], resultsTip[iHae, iC + 1], desCi)
            aeroshellTop = TrackCommon.Lerp(axisCi[iC], resultsRadiusTop[iHae, iC],
                                            axisCi[iC + 1], resultsRadiusTop[iHae, iC + 1], desCi)

        else:
            b1 = TrackCommon.Bilerp(axisCi[iC], axisCi[iC + 1], axisHae[iHae], axisHae[iHae + 1],
                                    resultsB1[iHae, iC], resultsB1[iHae + 1, iC],
                                    resultsB1[iHae + 1, iC + 1], resultsB1[iHae, iC + 1], desCi, desHae)
            b2 = TrackCommon.Bilerp(axisCi[iC], axisCi[iC + 1], axisHae[iHae], axisHae[iHae + 1],
                                    resultsB2[iHae, iC], resultsB2[iHae + 1, iC],
                                    resultsB2[iHae + 1, iC + 1], resultsB2[iHae, iC + 1], desCi, desHae)
            cTip = TrackCommon.Bilerp(axisCi[iC], axisCi[iC + 1], axisHae[iHae], axisHae[iHae + 1],
                                    resultsTip[iHae, iC], resultsTip[iHae + 1, iC],
                                    resultsTip[iHae + 1, iC + 1], resultsTip[iHae, iC + 1], desCi, desHae)
            aeroshellTop = TrackCommon.Bilerp(axisCi[iC], axisCi[iC + 1], axisHae[iHae], axisHae[iHae + 1],
                                    resultsRadiusTop[iHae, iC], resultsRadiusTop[iHae + 1, iC],
                                    resultsRadiusTop[iHae + 1, iC + 1], resultsRadiusTop[iHae, iC + 1], desCi, desHae)

        area = b1 * (cRoot + desCi) / 2 + b2 * (desCi + cTip) / 2
        print(TrackCommon.StringPad('Inboard half-span:       ', b1, 3, 8) + ' m')
        print(TrackCommon.StringPad('Outboard half-span:      ', b2, 3, 8) + ' m')
        print(TrackCommon.StringPad('Root chord:              ', cRoot, 3, 8) + ' m')
        print(TrackCommon.StringPad('Mid chord:               ', desCi, 3, 8) + ' m')
        print(TrackCommon.StringPad('Tip chord:               ', cTip, 3, 8) + ' m')
        print(TrackCommon.StringPad('Inboard taper ratio:     ', desCi / cRoot, 3, 8) + ' m')
        print(TrackCommon.StringPad('Outboard taper ratio:    ', cTip / desCi, 3, 8) + ' m')
        print(TrackCommon.StringPad('Aeroshell bottom radius: ', rAeroBottom, 3, 8) + ' m')
        print(TrackCommon.StringPad('Aeroshell top radius:    ', aeroshellTop, 3, 8) + ' m')
        print(TrackCommon.StringPad('Aeroshell height:        ', desHae, 3, 8) + ' m')
        print(TrackCommon.StringPad('Half area:               ', area, 3, 8) + ' m2')

def HighlightOperationalPoints(fileCl, fileCd, ReOperational, ClOperational, title=None):
    # Check the input for consistency
    if len(ReOperational) != len(ClOperational):
        raise ValueError("length of operational Reynolds numbers does not match " +
            "the number of operational Cl values")
        
    # Load the lookup tables
    lookupCl, lookupCd = TrackIO.LoadAerodynamicReynoldsData(fileCl, fileCd)
    axisAlpha, axisRe, mapCl = lookupCl.getPoints()
    _, _, mapCd = lookupCd.getPoints()
    
    # Finding the location of all operational points
    alphaOperational = np.zeros([len(ReOperational)])
    CdOperational = np.zeros(alphaOperational.shape)
    
    reverseCl = TrackLookup.Lookup2DReverse(*lookupCl.getPoints())
    labelOperational = []
    
    for i in range(0, len(ReOperational)):
        alphaOperational[i] = reverseCl.find(ClOperational[i], ReOperational[i])
        CdOperational[i] = lookupCd(alphaOperational[i], ReOperational[i])
        labelOperational.append('Cl = ' + str(round(ClOperational[i], 3)) + 
            '\nRe = ' + str(round(ReOperational[i] / 1e6, 3)) + 'M')
        
    ClOperational = np.asarray(ClOperational)
    CdOperational = np.asarray(CdOperational)
    ClOverCdOperational = ClOperational / CdOperational
        
    # Plotting Cl data
    cmap = plt.get_cmap('gnuplot2')
    
    fig = plt.figure()
    axCl = fig.add_subplot(131)
    
    for iRe in range(0, len(axisRe)):
        axCl.plot(axisAlpha, mapCl[:, iRe], color=cmap((iRe + 1) / len(axisRe)))
        axCl.scatter(alphaOperational, ClOperational)
        
        for iOp in range(0, len(ReOperational)):
            axCl.text(alphaOperational[iOp], ClOperational[iOp], labelOperational[iOp])
    
    axCl.grid(True)
    
    # Plotting Cd data
    axCd = fig.add_subplot(132, sharex=axCl)
    
    for iRe in range(0, len(axisRe)):
        axCd.plot(axisAlpha, mapCd[:, iRe], color=cmap((iRe + 1) / len(axisRe)))
        axCd.scatter(alphaOperational, CdOperational)
        
        for iOp in range(0, len(ReOperational)):
            axCd.text(alphaOperational[iOp], CdOperational[iOp], labelOperational[iOp])
    
    axCd.grid(True)
    
    # Plotting Cl/Cd data
    axClOverCd = fig.add_subplot(133, sharex=axCl)
    ClOverCd = mapCl / mapCd
    
    for iRe in range(0, len(axisRe)):
        axClOverCd.plot(axisAlpha, ClOverCd[:, iRe], color=cmap((iRe + 1) / len(axisRe)))
        axClOverCd.scatter(alphaOperational, ClOverCdOperational)
        
        for iOp in range(0, len(ReOperational)):
            axClOverCd.text(alphaOperational[iOp], ClOverCdOperational[iOp], labelOperational[iOp])
            
    axClOverCd.grid(True)
            
    if not (title is None):
        fig.suptitle(title)
    
#DetermineLengths(18.0, 4.58 / 2, 7.0, 4.0, 3.04332, 4.93848)
#HighlightOperationalPoints('./data/aerodynamicPerformance/v2/ClClimber.csv',
#                           './data/aerodynamicPerformance/v2/CdClimber.csv',
#                           [10**6.7, 10**7.1, 10**6.4, 10**6.9, 10**6.2, 10**6.7], 
#                           [0.5796, 0.5796, 0.3519, 0.3519, 0.3000, 0.3000],
#                           'Climber')
#HighlightOperationalPoints('./data/aerodynamicPerformance/v2/ClCruiser.csv',
#                           './data/aerodynamicPerformance/v2/CdCruiser.csv',
#                           [10**6.7, 10**7.1, 10**6.4, 10**6.9, 10**6.2, 10**6.7], 
#                           [0.5796, 0.5796, 0.3519, 0.3519, 0.3000, 0.3000],
#                           'Cruiser')
#HighlightOperationalPoints('./data/aerodynamicPerformance/v3/Cl.csv',
#                           './data/aerodynamicPerformance/v3/Cd.csv',
#                           [10**6.7, 10**7.1, 10**6.4, 10**6.9, 10**6.2, 10**6.7], 
#                           [0.5796, 0.5796, 0.3519, 0.3519, 0.3000, 0.3000],
#                           'OptCruiser')
HighlightOperationalPoints('./data/aerodynamicPerformance/v4/Cl.csv',
                           './data/aerodynamicPerformance/v4/Cd.csv',
                           [10**6.7, 10**7.1, 10**6.4, 10**6.9, 10**6.2, 10**6.7], 
                           [0.5796, 0.5796, 0.3519, 0.3519, 0.3000, 0.3000],
                           'OptCruiser2')
HighlightOperationalPoints('./data/aerodynamicPerformance/v4/oldCl.csv',
                           './data/aerodynamicPerformance/v4/oldCd.csv',
                           [10**6.7, 10**7.1, 10**6.4, 10**6.9, 10**6.2, 10**6.7], 
                           [0.5796, 0.5796, 0.3519, 0.3519, 0.3000, 0.3000],
                           'OptCruiserOld')