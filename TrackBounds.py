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
              cmap='gnuplot2', contours=None, forceNormMin=None, forceNormMax=None):
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

    axImage.imshow(data, extent=[xMin, xMax, yMin, yMax], norm=norm,
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

    return qInf, vInf, alpha, thrust, PReq, vInfOvera

def PlotMaps(axisHeight, axisDeltaV, maps):
    fig = plt.figure()
    bordersQInfData, bordersQInfLegend = ImageAxes(0.0, 0.5, 0.5, 1.0)
    bordersAlphaData, bordersAlphaLegend = ImageAxes(0.0, 0.5, 0.0, 0.5)
    bordersPReqData, bordersPReqLegend = ImageAxes(0.5, 1.0, 0.5, 1.0)
    bordersRawVInfData, bordersVInfLegend = ImageAxes(0.5, 1.0, 0.0, 0.5)

    axQInfData = fig.add_axes(bordersQInfData)
    axQInfLegend = fig.add_axes(bordersQInfLegend)

    PlotImage(fig, axQInfData, axQInfLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', np.log10(maps[0]), r'$log_10(q_\infty)\;[Pa]$',
        contours=[2.0, 2.5, 3.0, 3.5, 4.0, 4.5], forceNormMin=0)

    axAlphaData = fig.add_axes(bordersAlphaData)
    axAlphaLegend = fig.add_axes(bordersAlphaLegend)

    PlotImage(fig, axAlphaData, axAlphaLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', maps[2], r'$\alpha\;[\degree]$', cmap='gnuplot2_r',
        contours=[0, 0.5, 1.0, 2.5, 5.0, 7.5])

    axPReqData = fig.add_axes(bordersPReqData)
    axPReqLegend = fig.add_axes(bordersPReqLegend)

    PlotImage(fig, axPReqData, axPReqLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', maps[4] / 1e3, r'$P_\mathrm{req}\;[kW]$',
        cmap='gnuplot_r', contours=[0, 10, 20, 25, 30, 35, 40, 45, 50],
        forceNormMin=0, forceNormMax=100)

    axVInfData = fig.add_axes(bordersRawVInfData)
    axVInfLegend = fig.add_axes(bordersVInfLegend)

    PlotImage(fig, axVInfData, axVInfLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', maps[5], r'$V_\infty\/a;[-]$',
        contours=[0.2, 0.4, 0.5, 0.6, 0.7])

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
        ax.contourf(deltaVMesh, heightMesh / 1e3, maps[iSeverity][0], levels=[qInfMin, qInfMax], 
                    zdir='z', offset=axisSeverity[iSeverity], colors='r', alpha=0.2)
        ax.contourf(deltaVMesh, heightMesh / 1e3, maps[iSeverity][2], levels=[alphaMin, alphaMax],
                    zdir='z', offset=axisSeverity[iSeverity], colors='g', alpha=0.2)
        ax.contourf(deltaVMesh, heightMesh / 1e3, maps[iSeverity][4], levels=[PReqMin, PReqMax],
                    zdir='z', offset=axisSeverity[iSeverity], colors='g', alpha=0.2)
        
    ax.set_zlim(severityMin, severityMax)
    ax.set_xlabel(r'$\Delta V\;[m/s]$')
    ax.set_ylabel(r'$h\;[km]$')
    ax.set_zlabel(r'$C_\mathrm{severity}$')

def __testPlotMaps__(severity):
    atm = Atmosphere.Atmosphere()
    axisHeight = np.linspace(10, 70, 250) * 1e3
    axisDeltaV = np.linspace(-40, 40, 250)
    lookupCl, lookupCd = TrackCommon.LoadAerodynamicData("data/aerodynamicPerformance/Cl.csv",
                                                         "data/aerodynamicPerformance/Cd.csv")
    maps = GenerateMaps(axisHeight, axisDeltaV, 0, 0, severity, 700*8.8, 35, 0, 
                        lookupCl, lookupCd, atm)
    PlotMaps(axisHeight, axisDeltaV, maps)

def __testPlotMaps3D__():
    atm = Atmosphere.Atmosphere()
    lookupCl, lookupCd = TrackCommon.LoadAerodynamicData("data/aerodynamicPerformance/Cl.csv",
                                                         "data/aerodynamicPerformance/Cd.csv")
    PlotMaps3D(30e3, 75e3, 150, -50, 50, 150, -1.5, 1.5, 5, 0, 0, 700*8.8, 35, 0,
        lookupCl, lookupCd, atm, -8, 8, 1*10**2.5, 1*10**4.5, 0, 30e3, 0, 0.7)

#__testPlotMaps__(-1.0)
__testPlotMaps3D__()
