#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 12:22:07 2016

@author: MaxHenger
"""

import Atmosphere
import TrackBounds
import TrackClimbOptimize
import TrackCommon
import TrackContour
import TrackSettings

import numpy as np
import matplotlib.pyplot as plt

def PlotCruiseMaps(axisHeight, axisDeltaV, maps, title=None):
    fig = plt.figure()
    bordersQInfData, bordersQInfLegend = TrackCommon.ImageAxes(0.0, 0.333, 0.00, 1.0)
    bordersAlphaData, bordersAlphaLegend = TrackCommon.ImageAxes(0.333, 0.666, 0.0, 0.5)
    bordersRawVInfData, bordersVInfLegend = TrackCommon.ImageAxes(0.666, 1.0, 0.0, 0.5)
    bordersPReqData, bordersPReqLegend = TrackCommon.ImageAxes(0.333, 0.666, 0.5, 1.0)
    bordersReData, bordersReLegend = TrackCommon.ImageAxes(0.666, 1.0, 0.5, 1.0)

    axQInfData = fig.add_axes(bordersQInfData)
    axQInfLegend = fig.add_axes(bordersQInfLegend)

    TrackCommon.PlotImage(fig, axQInfData, axQInfLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', np.log10(maps["qInf"]), r'$log_10(q_\infty)\;[Pa]$',
        contours=[2.0, 2.5, 3.0, 3.5, 4.0, 4.5], forceNormMin=0)
    contourQInf = TrackContour.Contour()
    contourQInf.setData(axisDeltaV, axisHeight, np.log10(maps["qInf"]), 2.0, 4.0)

    for i in range(0, contourQInf.getNumContours()):
        curContour = contourQInf.getContour(i)
        curIncluder = curContour.getIncluder()

        axQInfData.plot(curIncluder[:, 0], curIncluder[:, 1] / 1e3, 'g')

        for j in range(0, curContour.getNumExcluders()):
            curExcluder = curContour.getExcluder(j)

            axQInfData.plot(curExcluder[:, 0], curExcluder[:, 1] / 1e3, 'r')

    axAlphaData = fig.add_axes(bordersAlphaData)
    axAlphaLegend = fig.add_axes(bordersAlphaLegend)

    TrackCommon.PlotImage(fig, axAlphaData, axAlphaLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', maps["alpha"], r'$\alpha\;[\degree]$', cmap='gnuplot2_r',
        contours=[0, 0.5, 1.0, 2.5, 5.0, 7.5])

    axPReqData = fig.add_axes(bordersPReqData)
    axPReqLegend = fig.add_axes(bordersPReqLegend)

    TrackCommon.PlotImage(fig, axPReqData, axPReqLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', maps["PReq"] / 1e3, r'$P_\mathrm{req}\;[kW]$',
        cmap='gnuplot_r', contours=[0, 10, 20, 25, 30, 35, 40, 45, 50],
        forceNormMin=0, forceNormMax=100)

    axVInfData = fig.add_axes(bordersRawVInfData)
    axVInfLegend = fig.add_axes(bordersVInfLegend)

    TrackCommon.PlotImage(fig, axVInfData, axVInfLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', maps["vInfOvera"], r'$V_\infty\/a;[-]$',
        contours=[0.2, 0.4, 0.5, 0.6, 0.7])

    axReData = fig.add_axes(bordersReData)
    axReLegend = fig.add_axes(bordersReLegend)

    TrackCommon.PlotImage(fig, axReData, axReLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', np.log10(maps["Re"]),
        r'$\mathrm{log}_{10} ( Re ) \; [-]$', forceNormMin=0)

    if title != None:
        fig.suptitle(title)

def PlotReportCruiseMap(axisHeight, axisDeltaV, maps, includeLegend=False):
    fig = plt.figure()
    boundsData, boundsLegend = TrackCommon.ImageAxes(0.07, 0.93, 0.08, 1.0)
    axData = fig.add_axes(boundsData)
    axLegend = fig.add_axes(boundsLegend)
    
    TrackCommon.PlotImage(fig, axData, axLegend, axisDeltaV, r'$V_{I,\mathrm{hor}} \; [m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', maps["PReq"] / 1e3, r'$P_A \; [kW]$',
        cmap='gnuplot_r', forceNormMin=0, forceNormMax=100, fontsize=15)
    
    cmap = plt.get_cmap('jet')
    lines = []
    names = []

    PReqList = [20e3, 26e3, 32e3, 36e3, 40e3]
    for iPReq in range(len(PReqList)):
        PReq = PReqList[iPReq]
        contour = TrackContour.Contour()
        contour.combineData(axisDeltaV, axisHeight, [maps['PReq'],
            maps['qInf'], maps['vInfOvera'], maps['alpha']],
            [0, 200, 0, -8.0], [PReq, 1e10, 0.65, 8.0])
        
        for i in range(contour.getNumContours()):
            curContour = contour.getContour(i)
            includer = curContour.getIncluder()
            curColor = cmap((iPReq + 1) / len(PReqList))
            
            newLine, = axData.plot(includer[:,0], includer[:,1] / 1e3,
                                   color=curColor)
            
            if i == 0:
                lines.append(newLine)
                names.append(r'$P_A = ' + str(round(PReq / 1e3, 1)) + ' kW$')
                
            for j in range(curContour.getNumExcluders()):
                excluder = curContour.getExcluder(j)
                axData.plot(excluder[:, 0], excluder[:, 1] / 1e3,
                            linestyle='--', color=curColor)
       
    if includeLegend:
        axData.legend(lines, names, loc=4)
    
def PlotCombinedMaps(heightMin, heightMax, heightNum, deltaVMin, deltaVMax,
                     deltaVNum, PRequired, severity):
    # Retrieve the climbing maps. The only real important variables are the
    # power required and the horizontal velocity at which it occurs. So only
    # those two things are plotted.

    # Generate axes and often-used variables
    axisHeight = np.linspace(heightMin, heightMax, heightNum)
    axisDeltaV = np.linspace(deltaVMin, deltaVMax, deltaVNum)
    settings = TrackSettings.Settings()
    atmosphere = Atmosphere.Atmosphere()

    if not TrackCommon.IsArray(PRequired):
        PRequired = [PRequired]

    # Generate the cruise maps
    cruiseDict = TrackBounds.GenerateCruiseMaps(axisHeight, axisDeltaV,
        settings.latitude, settings.longitude, severity, settings.W, settings.S,
        settings.inclination, settings.lookupCl, settings.lookupCd, atmosphere,
        settings.reynoldsLength)

    # Create all required contours
    contours = []

    for i in range(0, len(PRequired)):
        # - Generate the climbing maps
        climbDict = TrackClimbOptimize.GenerateAscentMaps(axisHeight, axisDeltaV,
            settings.W, settings.S, settings.inclination, settings.lookupCl,
            settings.lookupCd, atmosphere, settings.qInfMin, settings.qInfMax,
            settings.alphaMin, settings.alphaMax, 0, PRequired[i],
            settings.latitude, settings.longitude, settings.reynoldsLength,
            severity, storeResults=False)

        newContour = TrackContour.Contour()
        newContour.combineData(axisDeltaV, axisHeight,
            [cruiseDict['qInf'], climbDict['qInf'], cruiseDict['PReq'],
             climbDict['PReq'], cruiseDict['alpha'], climbDict['alpha']],
            [settings.qInfMin, settings.qInfMin, 0,
             0, settings.alphaMin, settings.alphaMin],
            [settings.qInfMax, settings.qInfMax, PRequired[i],
             PRequired[i], settings.alphaMax, settings.alphaMax])

        contours.append(newContour)

    fig = plt.figure()
    bordersMap, bordersLegend = TrackCommon.ImageAxes(0.1, 0.9, 0.1, 0.9)
    axMap = fig.add_axes(bordersMap)
    axLegend = fig.add_axes(bordersLegend)
    TrackCommon.PlotImage(fig, axMap, axLegend, axisDeltaV, r'$\Delta V \; [m/s]$',
        axisHeight / 1e3, r'$h \; [km]$', cruiseDict['PReq'] / 1e3,
        r'$P_\mathrm{req} \; [kW]$')

    contourColormap = plt.get_cmap('GnBu')
    contourLines = []
    contourNames = []

    for i in range(0, len(contours)):
        curContour = contours[i]

        for j in range(0, curContour.getNumContours()):
            curSubContour = curContour.getContour(j)
            includer = curSubContour.getIncluder()
            contourLine, = axMap.plot(includer[:, 0], includer[:, 1] / 1e3,
                color=contourColormap((i + 1) / len(contours)))

            if j == 0:
                contourLines.append(contourLine)

            for k in range(0, curSubContour.getNumExcluders()):
                excluder = curSubContour.getExcluder(k)

                axMap.plot(excluder[:, 0], excluder[:, 1] / 1e3,
                    color=contourColormap((i + 1) / len(contours)),
                    linestyle='--')

        contourNames.append(r'$P_\mathrm{req} \; = \; ' +
                            str(round(PRequired[i] / 1e3, 0)) + ' kW$')

    axMap.legend(contourLines, contourNames)

    # Go through the combined map and see where valid regions of cruise exist


def PlotCruiseMaps3D(heightMin, heightMax, heightNum, deltaVMin, deltaVMax,
        deltaVNum, severityMin, severityMax, severityNum, latitude, longitude,
        W, S, inclination, lookupCl, lookupCd, atm, alphaMin, alphaMax, qInfMin,
        qInfMax, PReqMin, PReqMax, VInfOveraMin, VInfOveraMax):
    axisHeight = np.linspace(heightMin, heightMax, heightNum)
    axisDeltaV = np.linspace(deltaVMin, deltaVMax, deltaVNum)
    axisSeverity = np.linspace(severityMin, severityMax, severityNum)

    maps = []
    ascentMaps = []

    for severity in axisSeverity:
        maps.append(TrackBounds.GenerateCruiseMaps(axisHeight, axisDeltaV,
            latitude, longitude, severity, W, S, inclination, lookupCl,
            lookupCd, atm, settings.reynoldsLength))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #curMaps = maps[0]
    #PlotMaps(axisHeight, axisDeltaV, curMaps)
    deltaVMesh, heightMesh = np.meshgrid(axisDeltaV, axisHeight)

    for iSeverity in range(0, severityNum):
        ax.contourf(deltaVMesh, heightMesh / 1e3, maps[iSeverity]['qInf'], levels=[qInfMin, qInfMax],
                    zdir='z', offset=axisSeverity[iSeverity], colors='r', alpha=0.2)
        ax.contourf(deltaVMesh, heightMesh / 1e3, maps[iSeverity]['alpha'], levels=[alphaMin, alphaMax],
                    zdir='z', offset=axisSeverity[iSeverity], colors='g', alpha=0.2)
        ax.contourf(deltaVMesh, heightMesh / 1e3, maps[iSeverity]['PReq'], levels=[PReqMin, PReqMax],
                    zdir='z', offset=axisSeverity[iSeverity], colors='g', alpha=0.2)

    ax.set_zlim(severityMin, severityMax)
    ax.set_xlabel(r'$\Delta V\;[m/s]$')
    ax.set_ylabel(r'$h\;[km]$')
    ax.set_zlabel(r'$C_\mathrm{severity}$')

def PlotCruiseContour(ax, contour, z, colorIncluder='g', linestyleIncluder='-',
                      colorExcluder='r'):
    for i in range(0, contour.getNumContours()):
        curContour = contour.getContour(i)
        curIncluder = curContour.getIncluder()

        ax.plot(curIncluder[:, 0], curIncluder[:, 1],
            np.repeat(z, len(curIncluder)), color=colorIncluder,
            linestyle=linestyleIncluder)

def PlotCruiseContour3D(heightMin, heightMax, heightNum, deltaVMin, deltaVMax,
        deltaVNum, severityMin, severityMax, severityNum, latitude, longitude,
        W, S, inclination, lookupCl, lookupCd, atm, alphaMin, alphaMax, qInfMin,
        qInfMax, PReqMin, PReqMax, vInfOveraMin, vInfOveraMax, includeAll=False):
    axisHeight = np.linspace(heightMin, heightMax, heightNum)
    axisDeltaV = np.linspace(deltaVMin, deltaVMax, deltaVNum)
    axisSeverity = np.linspace(severityMin, severityMax, severityNum)

    maps = []

    for severity in axisSeverity:
        maps.append(TrackBounds.GenerateCruiseMaps(axisHeight, axisDeltaV,
            latitude, longitude, severity, W, S, inclination, lookupCl,
            lookupCd, atm, settings.reynoldsLength))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    cmap = plt.get_cmap('jet')

    for i in range(0, len(axisSeverity)):
        if includeAll:
            contourAlpha = TrackContour.Contour()
            contourAlpha.setData(axisDeltaV, axisHeight / 1e3, maps[i]['alpha'], alphaMin, alphaMax)
            PlotCruiseContour(ax, contourAlpha, axisSeverity[i], colorIncluder='r', linestyleIncluder='--')

            contourQInf = TrackContour.Contour()
            contourQInf.setData(axisDeltaV, axisHeight / 1e3, maps[i]['qInf'], qInfMin, qInfMax)
            PlotCruiseContour(ax, contourQInf, axisSeverity[i], colorIncluder='g', linestyleIncluder='--')

            contourPReq = TrackContour.Contour()
            contourPReq.setData(axisDeltaV, axisHeight / 1e3, maps[i]['PReq'], PReqMin, PReqMax)
            PlotCruiseContour(ax, contourPReq, axisSeverity[i], colorIncluder='b', linestyleIncluder='--')

            contourVInfOvera = TrackContour.Contour()
            contourVInfOvera.setData(axisDeltaV, axisHeight / 1e3, maps[i]['vInfOvera'], vInfOveraMin, vInfOveraMax)
            PlotCruiseContour(ax, contourVInfOvera, axisSeverity[i], colorIncluder='k', linestyleIncluder='--')

        contourCombined = TrackContour.Contour()
        contourCombined.combineData(axisDeltaV, axisHeight / 1e3,
            [maps[i]['alpha'], maps[i]['qInf'], maps[i]['PReq'], maps[i]['vInfOvera']],
            [alphaMin, qInfMin, PReqMin, vInfOveraMin],
            [alphaMax, qInfMax, PReqMax, vInfOveraMax])

        PlotCruiseContour(ax, contourCombined, axisSeverity[i], colorIncluder=cmap((i + 1) / len(axisSeverity)))

    ax.set_xlabel(r'$\Delta V\;[m/s]$')
    ax.set_ylabel(r'$h\;[km]$')
    ax.set_zlabel(r'$C_\mathrm{severity}$')

def __testPlotMaps__(severity):
    atm = Atmosphere.Atmosphere()
    axisHeight = np.linspace(10, 75, 125) * 1e3
    axisDeltaV = np.linspace(-80, 50, 125)
    settings = TrackSettings.Settings()
    maps = TrackBounds.GenerateCruiseMaps(axisHeight, axisDeltaV,
        settings.latitude, settings.longitude, severity, settings.W, settings.S,
        settings.inclination, settings.lookupCl, settings.lookupCd, atm,
        settings.reynoldsLength)
    PlotCruiseMaps(axisHeight, axisDeltaV, maps, 'severity = ' + str(round(severity, 3)))

def __testPlotReportMap__(severity, minDeltaV, maxDeltaV, numDeltaV, includeLegend=False):
    atm = Atmosphere.Atmosphere()
    axisHeight = np.linspace(30, 80, 125) * 1e3
    axisDeltaV = np.linspace(minDeltaV, maxDeltaV, numDeltaV)
    settings = TrackSettings.Settings()
    maps = TrackBounds.GenerateCruiseMaps(axisHeight, axisDeltaV,
        settings.latitude, settings.longitude, severity, settings.W, settings.S,
        settings.inclination, settings.lookupCl, settings.lookupCd, atm,
        settings.reynoldsLength)
    PlotReportCruiseMap(axisHeight, axisDeltaV, maps, includeLegend)
    
def __testPlotMaps3D__():
    atm = Atmosphere.Atmosphere()
    settings = TrackSettings.Settings()
    PlotCruiseMaps3D(30e3, 80e3, 75, -50, 50, 75, -1.5, 1.5, 4, settings.latitude,
        settings.longitude, settings.W, settings.S, settings.inclination,
        settings.lookupCl, settings.lookupCd, atm, -8, 8, 200, 1*10**10, 0, 32e3,
        0, 0.7)

def __testPlotContour3D__():
    atm = Atmosphere.Atmosphere()
    settings = TrackSettings.Settings()
    PlotCruiseContour3D(30e3, 80e3, 75, -50, 50, 75, -2.5, 2.5, 11,
        settings.latitude, settings.longitude, settings.W, settings.S,
        settings.inclination, settings.lookupCl, settings.lookupCd, atm, -8, 8,
        200, 1*10**10, 0, 32e3, 0, 0.7, includeAll=False)

def __testCombinedMaps__(severity):
    PlotCombinedMaps(30e3, 80e3, 60, -80, 40, 70, [20e3, 30e3, 40e3], severity)

#__testPlotMaps__(1.6)
#__testPlotReportMap__(-1.55, -20, 120, 125)
#__testPlotReportMap__(0.0, -40, 80, 125)
#__testPlotReportMap__(1.8, -100, 40, 125, includeLegend=True)
#__testPlotMaps__(0.0)
#__testPlotMaps__(1.6)
#__testPlotMaps3D__()
#__testPlotContour3D__()
#__testCombinedMaps__(-1.6)
#__testCombinedMaps__(0.0)
#__testCombinedMaps__(1.6)
