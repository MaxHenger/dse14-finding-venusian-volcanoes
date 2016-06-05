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

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def ImageAxes(xl, xr, yb, yt, colorbarSeperation=0.05, colorbarWidth=0.03, border=0.1):
    return [ # The plot's axes
        xl + border / 2,
        yb + border / 2,
        xr - xl - colorbarWidth - colorbarSeperation - border,
        yt - yb - border
    ], \
    [ # The colorbar's axes
        xr - colorbarWidth - colorbarSeperation / 2 - border / 2,
        yb + border / 2,
        colorbarWidth,
        yt - yb - border
    ]

def PlotImage(fig, axImage, axColorbar, xAxis, xLabel, yAxis, yLabel, data, dataLabel,
              cmap='gnuplot2', contours=None, forceNormMin=None, forceNormMax=None,
              extentAdjust=False):
    # Determine extent of axes
    xMin = min(xAxis)
    xMax = max(xAxis)
    yMin = min(yAxis)
    yMax = max(yAxis)

    # Determine normalization
    vMin = 0
    vMax = 0

    if forceNormMin == None:
        vMin = np.min(data)
    else:
        vMin = forceNormMin

    if forceNormMax == None:
        vMax = np.max(data)
    else:
        vMax = forceNormMax

    norm = mpl.colors.Normalize(vMin, vMax)

    # Determine colormap
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)

    # Determine extent. shift by half a pixel leftward and downward if required.
    # This will assume all data is equally spaced
    extent = [xMin, xMax, yMin, yMax]

    if extentAdjust == True:
        deltaX = (xAxis[1] - xAxis[0]) / 2.0
        deltaY = (yAxis[1] - yAxis[0]) / 2.0
        xMin -= deltaX
        xMax += deltaX
        yMin -= deltaY
        yMax += deltaY

    axImage.imshow(data, extent=extent, norm=norm,
                   aspect='auto', origin='lower', cmap=cmap)
    axImage.set_xlabel(xLabel)
    axImage.set_ylabel(yLabel)
    axImage.grid(True)

    if contours != None:
        ct = axImage.contour(xAxis, yAxis, data, contours, colors='k')
        axImage.clabel(ct)

    # Add the colorbar
    cbb = mpl.colorbar.ColorbarBase(axColorbar, cmap=cmap, norm=norm)
    cbb.set_label(dataLabel, rotation=90, fontsize=14)

def GenerateMaps(axisHeight, axisDeltaV, latitude, longitude, relativeSeverity,
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

    velocityZonal = atm.velocityZonal(axisHeight, latitude, longitude)
    density = atm.density(axisHeight, latitude, longitude)
    speedOfSound = atm.speedOfSound(axisHeight, latitude, longitude)

    # Adjust the velocity map with respect to the 'relative severity'
    if relativeSeverity >= 0.0:
        velocityZonal = velocityZonal[1] + (velocityZonal[2] - velocityZonal[1]) * relativeSeverity
        density = density[1] + (density[2] - density[1]) * relativeSeverity
    else:
        velocityZonal = velocityZonal[1] + (velocityZonal[1] - velocityZonal[0]) * relativeSeverity
        density = density[1] + (density[1] - density[0]) * relativeSeverity

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

def PlotMaps(axisHeight, axisDeltaV, maps, title=None):
    fig = plt.figure()
    bordersQInfData, bordersQInfLegend = ImageAxes(0.0, 0.5, 0.5, 1.0)
    bordersAlphaData, bordersAlphaLegend = ImageAxes(0.0, 0.5, 0.0, 0.5)
    bordersPReqData, bordersPReqLegend = ImageAxes(0.5, 1.0, 0.5, 1.0)
    bordersRawVInfData, bordersVInfLegend = ImageAxes(0.5, 1.0, 0.0, 0.5)

    axQInfData = fig.add_axes(bordersQInfData)
    axQInfLegend = fig.add_axes(bordersQInfLegend)

    PlotImage(fig, axQInfData, axQInfLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
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

    PlotImage(fig, axAlphaData, axAlphaLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', maps["alpha"], r'$\alpha\;[\degree]$', cmap='gnuplot2_r',
        contours=[0, 0.5, 1.0, 2.5, 5.0, 7.5])

    axPReqData = fig.add_axes(bordersPReqData)
    axPReqLegend = fig.add_axes(bordersPReqLegend)

    PlotImage(fig, axPReqData, axPReqLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', maps["PReq"] / 1e3, r'$P_\mathrm{req}\;[kW]$',
        cmap='gnuplot_r', contours=[0, 10, 20, 25, 30, 35, 40, 45, 50],
        forceNormMin=0, forceNormMax=100)

    axVInfData = fig.add_axes(bordersRawVInfData)
    axVInfLegend = fig.add_axes(bordersVInfLegend)

    PlotImage(fig, axVInfData, axVInfLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', maps["vInfOvera"], r'$V_\infty\/a;[-]$',
        contours=[0.2, 0.4, 0.5, 0.6, 0.7])

    if title != None:
        fig.suptitle(title)

def PlotMaps3D(heightMin, heightMax, heightNum, deltaVMin, deltaVMax, deltaVNum,
        severityMin, severityMax, severityNum, latitude, longitude, W, S,
        inclination, lookupCl, lookupCd, atm, alphaMin, alphaMax, qInfMin,
        qInfMax, PReqMin, PReqMax, VInfOveraMin, VInfOveraMax):
    axisHeight = np.linspace(heightMin, heightMax, heightNum)
    axisDeltaV = np.linspace(deltaVMin, deltaVMax, deltaVNum)
    axisSeverity = np.linspace(severityMin, severityMax, severityNum)

    maps = []

    for severity in axisSeverity:
        maps.append(GenerateMaps(axisHeight, axisDeltaV, latitude, longitude,
            severity, W, S, inclination, lookupCl, lookupCd, atm))

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

def PlotContour(ax, contour, z, colorIncluder='g', colorExcluder='r'):
    for i in range(0, contour.getNumContours()):
        curContour = contour.getContour(i)
        curIncluder = curContour.getIncluder()

        ax.plot(curIncluder[:, 0], curIncluder[:, 1],
            np.repeat(z, len(curIncluder)), color=colorIncluder)

        for j in range(0, curContour.getNumExcluders()):
            curExcluder = curContour.getExcluder(j)

            ax.plot(curExcluder[:, 0], curExcluder[:, 1],
                np.repeat(z, len(curExcluder)), color=colorExcluder)

def PlotContour3D(heightMin, heightMax, heightNum, deltaVMin, deltaVMax, deltaVNum,
        severityMin, severityMax, severityNum, latitude, longitude, W, S,
        inclination, lookupCl, lookupCd, atm, alphaMin, alphaMax, qInfMin,
        qInfMax, PReqMin, PReqMax, vInfOveraMin, vInfOveraMax, includeAll=False):
    axisHeight = np.linspace(heightMin, heightMax, heightNum)
    axisDeltaV = np.linspace(deltaVMin, deltaVMax, deltaVNum)
    axisSeverity = np.linspace(severityMin, severityMax, severityNum)

    maps = []

    for severity in axisSeverity:
        maps.append(GenerateMaps(axisHeight, axisDeltaV, latitude, longitude,
            severity, W, S, inclination, lookupCl, lookupCd, atm))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range(0, len(axisSeverity)):
        if includeAll:
            contourAlpha = TrackContour.Contour()
            contourAlpha.setData(axisDeltaV, axisHeight, maps[i]['alpha'], alphaMin, alphaMax)
            PlotContour(ax, contourAlpha, axisSeverity[i])
    
            contourQInf = TrackContour.Contour()
            contourQInf.setData(axisDeltaV, axisHeight, maps[i]['qInf'], qInfMin, qInfMax)
            PlotContour(ax, contourQInf, axisSeverity[i])
            
            contourPReq = TrackContour.Contour()
            contourPReq.setData(axisDeltaV, axisHeight, maps[i]['PReq'], PReqMin, PReqMax)
            PlotContour(ax, contourPReq, axisSeverity[i])
            
            contourVInfOvera = TrackContour.Contour()
            contourVInfOvera.setData(axisDeltaV, axisHeight, maps[i]['vInfOvera'], vInfOveraMin, vInfOveraMax)
            PlotContour(ax, contourVInfOvera, axisSeverity[i])
        
        contourCombined = TrackContour.Contour()
        contourCombined.combineData(axisDeltaV, axisHeight,
            [maps[i]['alpha'], maps[i]['qInf'], maps[i]['PReq'], maps[i]['vInfOvera']],
            [alphaMin, qInfMin, PReqMin, vInfOveraMin],
            [alphaMax, qInfMax, PReqMax, vInfOveraMax])
        
        PlotContour(ax, contourCombined, axisSeverity[i], colorIncluder='b')

def __testPlotMaps__(severity):
    atm = Atmosphere.Atmosphere()
    axisHeight = np.linspace(10, 75, 125) * 1e3
    axisDeltaV = np.linspace(-50, 50, 125)
    lookupCl, lookupCd = TrackCommon.LoadAerodynamicData("data/aerodynamicPerformance/Cl.csv",
                                                         "data/aerodynamicPerformance/Cd.csv")
    maps = GenerateMaps(axisHeight, axisDeltaV, 0, 0, severity, 700*8.8, 35, 0,
                        lookupCl, lookupCd, atm)
    PlotMaps(axisHeight, axisDeltaV, maps, 'severity = ' + str(round(severity, 3)))

def __testPlotMaps3D__():
    atm = Atmosphere.Atmosphere()
    lookupCl, lookupCd = TrackCommon.LoadAerodynamicData("data/aerodynamicPerformance/Cl.csv",
                                                         "data/aerodynamicPerformance/Cd.csv")
    PlotMaps3D(30e3, 80e3, 75, -50, 50, 75, -1.5, 1.5, 4, 0, 0, 700*8.8, 35, 0,
        lookupCl, lookupCd, atm, -8, 8, 1*10**2.5, 1*10**4.5, 0, 32e3, 0, 0.7)

def __testPlotContour3D__():
    atm = Atmosphere.Atmosphere()
    lookupCl, lookupCd = TrackCommon.LoadAerodynamicData("data/aerodynamicPerformance/Cl.csv",
                                                         "data/aerodynamicPerformance/Cd.csv")
    PlotContour3D(30e3, 80e3, 50, -50, 50, 50, -2.5, 2.5, 15, 0, 0, 700*8.8, 35,
        0, lookupCl, lookupCd, atm, -8, 8, 1*10**2.5, 1*10**4.5, 0, 35e3, 0, 0.7)

#__testPlotMaps__(-1.5)
#__testPlotMaps__(0.0)
#__testPlotMaps__(1.5)
#__testPlotMaps3D__()
__testPlotContour3D__()
