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

    return {"qInf": qInf, "vInf": vInf, "alpha": alpha,
            "thrust": thrust, "PReq": PReq, "vInfOvera": vInfOvera}

# __determineSingleContourHorizontal__ is a function that will look for map
# contour points along a single horizontal line. This function will not check
# the input for validity to decrease computational cost.
#
# Input:
#   - height: A single height value in the desired units of the user
#   - axisDeltaV: The deltaV axis to look through in the desired units
#   - mapHoriz: A horizontal slice (along the deltaV axis) to look through to
#       find the contours
#   - vMin: The minimum bound of the contour
#   - vMax: The maximum bound of the contour
#
# Output:
#   - pointsLower: The points (an array of [iDeltaV, deltaV, iHeight, height]
#       indices and coordinates) that constitute the lower bounds of the contour
#   - pointsUpper: The points (an arary of [iDeltaV, deltaV, iHeight, height]
#       coordinates) that constitute the upper bounds of the contour
def __determineSingleContourHorizontal__(iHeight, height, axisDeltaV, mapHor, vMin, vMax):
    pointsLower = []
    pointsUpper = []

    for i in range(0, len(axisDeltaV) - 1):
        if mapHor[i + 1] >= mapHor[i]:
            # Next point is larger than current point, check for an intersection
            # with the minimum and maximum values going from low to high.
            if mapHor[i] <= vMin:
                if mapHor[i + 1] > vMin:
                    # Found a transition from low across vMin to high
                    pointsLower.append([
                        i,
                        TrackCommon.Lerp(mapHor[i], axisDeltaV[i],
                                         mapHor[i + 1], axisDeltaV[i + 1], vMin),
                        iHeight,
                        height
                    ])
            elif mapHor[i] <= vMax:
                if mapHor[i + 1] > vMax:
                    # Found a transition from low across vMax to high
                    pointsUpper.append([
                        i,
                        TrackCommon.Lerp(mapHor[i], axisDeltaV[i],
                                         mapHor[i + 1], axisDeltaV[i + 1], vMax),
                        iHeight,
                        height
                    ])
        else:
            # Next point is smaller than or equal to the current point, check
            # for an intersection with the minimum and maximum values from high
            # to low.
            if mapHor[i] >= vMax:
                if mapHor[i + 1] < vMax:
                    # Found a transition from high across vMax to low
                    pointsUpper.append([
                        i,
                        TrackCommon.Lerp(mapHor[i], axisDeltaV[i],
                                         mapHor[i + 1], axisDeltaV[i + 1], vMax),
                        iHeight,
                        height
                    ])
            elif mapHor[i] >= vMin:
                if mapHor[i + 1] < vMin:
                    # Found a transition from high across vMin to low
                    pointsLower.append([
                        i,
                        TrackCommon.Lerp(mapHor[i], axisDeltaV[i],
                                         mapHor[i + 1], axisDeltaV[i + 1], vMin),
                        iHeight,
                        height
                    ])

    return pointsLower, pointsUpper

# __determineSingleContourVertical__ is a function that will look for map
# contour points along a single vertical line. This function will not check the
# input for validity to decrease computational cost
#
# Input:
#   - iHeightLow: The index of the lower row of map values
#   - heightLow: The height value corresponding to the lower row of map values
#   - iHeightUp: The index of the upper row of map values
#   - heightUp: The height value corresponding to the upper row of map values
#   - axisDeltaV: The deltaV axis to look through
#   - mapLow: The lower row of map values to look through for intersections
#   - mapUp: The upper row of map values to look through for intersections
#   - vMin: The minimum bound of the contour
#   - vMax: The maximum bound of the contour
#
# Output:
#   - pointsLower: The points (an array of [deltaV, height] coordinates) that
#       constitute the lower bounds of the contour
#   - pointsUpper: The points (an array of [deltaV, height] coordinates) that
#       constitute the upper bounds of the contour
def __determineSingleContourVertical__(iHeightLow, heightLow, iHeightUp, heightUp,
        axisDeltaV, mapLow, mapUp, vMin, vMax):
    pointsLower = []
    pointsUpper = []

    for i in range(0, len(axisDeltaV)):
        if mapUp[i] > mapLow[i]:
            # Next point is larger than the current point, check for an
            # intersection with the minimum and maximum values going from low to
            # high
            if mapLow[i] <= vMin:
                if mapUp[i] > vMax:
                    # Found a transition from low across vMin to high
                    pointsLower.append([
                        i,
                        axisDeltaV[i],
                        iHeightLow,
                        TrackCommon.Lerp(mapLow[i], heightLow,
                                         mapUp[i], heightUp, vMin)
                    ])
            elif mapLow[i] <= vMax:
                if mapUp[i] > vMax:
                    # Found a transition from low across vMax to high
                    pointsUpper.append([
                        i,
                        axisDeltaV[i],
                        iHeightLow,
                        TrackCommon.Lerp(mapLow[i], heightLow,
                                         mapUp[i], heightUp, vMax)
                    ])
        else:
            # Next point is smaller than the current point, check for an
            # intersection with the minimum and maximum values going from high
            # to low
            if mapLow[i] >= vMax:
                if mapUp[i] < vMax:
                    # Found a transition from high across vMax to low
                    pointsUpper.append([
                        i,
                        axisDeltaV[i],
                        iHeightLow,
                        TrackCommon.Lerp(mapLow[i], heightLow,
                                         mapUp[i], heightUp, vMax)
                    ])
            elif mapVer[i] >= vMin:
                if mapVer[i + 1] < vMin:
                    # Found a transition from high across vMin to low
                    pointsLower.append([
                        i,
                        axisDeltaV[i],
                        iHeightLow,
                        TrackCommon.Lerp(mapLow[i], heightLow,
                                         mapUp[i], heightUp, vMin)
                    ])

    return pointsLower, pointsUpper

def __determineSingleContour__(axisHeight, axisDeltaV, data, vMin, vMax):
    # Make sure the input data is valid
    if len(map.shape) != 2:
        raise ValueError("Expected the 'map' to be a 2D dataset")

    if map.shape[0] != len(axisHeight):
        raise ValueError("Expected the first 'map' axis to equal 'axisHeight' length")

    if map.shape[1] != len(axisDeltaV):
        raise ValueError("Expected the second 'map' axis to equal 'axisDeltaV' length")

    if vMin >= vMax:
        raise ValueError("Expected 'vMin' to be smaller than 'vMax'")

    # Start calculating the bounds of the contours. Start at the lower height
    # and slowly work upwards. First start with the initial set of lines to ease
    # the iteration process
    linesLower, linesUpper = __determineSingleContourHorizontal__(0,
        axisHeight[0], axisDeltaV, data[0, :], vMin, vMax)

    extents=[]

    if len(linesLower) != 0 or len(linesUpper) != 0:
        # Create a new extent entry
        newExtent = [0, []]

        if len(linesLower) == 0:
            newExtent[0] = linesUpper[0][3]
        else:
            newExtent[0] = linesLower[0][3]

        if len(linesLower) != len(linesUpper):
            if abs(len(linesLower) - len(linesUpper)) > 1:
                raise RuntimeError("Expected length of 'linesLower' to differ " +
                    "by a maximum of one from 'linesUpper'")

            if len(linesLower) < len(linesUpper):
                # There is one upper-section at the left-most side
                newExtent[1].append([axisDeltaV[0], linesUpper[0][1]])

                for i in range(0, len(linesLower)):
                    newExtent[1].append([linesLower[i][1], linesUpper[i + 1][1]])
            else:
                # There is one lower-section to the right-most side
                for i in range(0, len(linesUpper)):
                    newExtent[1].append([linesLower[i][1], linesUpper[i][1]])

                newExtent[1].append([linesLower[-1][1], axisDeltaV[-1]])
        else:
            if linesUpper[0][1] > linesLower[0][1]:
                for i in range(0, len(linesLower)):
                    newExtent[1].append([linesLower[i][1], linesUpper[i][1]])
            else:
                # There is one upper-section to the left-most side and one
                # lower-section to the rightmost side
                newExtent[1].append([axisDeltaV[0], linesUpper[0][1]])

                for i in range(1, len(linesLower)):
                    newExtent[1].append([linesLower[i - 1][1], linesUpper[i][1]])

                newExtent[1].append([linesLower[-1][1], axisDeltaV[-1]])

        extents.append(newExtent)

    # Turn the initial coordinates into line arrays
    linesLower = [linesLower]
    linesUpper = [linesUpper]

    # Work through all the other contour lines
    for iHeight in range(1, len(axisHeight)):
        newLowerVer, newUpperVer = __determineSingleContourVertical__(i - 1,
            axisHeight[i - 1], i, axisHeight[i], axisDeltaV, data[i - 1, :],
            data[i, :], vMin, vMax)

        newLowerHor, newUpperHor = __determineSingleContourHorizontal__(i,
            axisHeight[i], axisDeltaV, data[i, :], vMin, vMax)

        # Assemble the slices for more compact code
        height = [(axisHeight[i - 1] + axisHeight[i]) / 2.0, axisHeight[i]]
        newLowerCombined = [newLowerVer, newLowerHor]
        newUpperCombined = [newUpperVer, newUpperHor]

        # See how the new lower vertical contour boundaries fit to the existing
        # lines. The current algorithm simply checks if the indices are close
        # enough
        for iSection in range(0, 2):
            # TODO: In case of incorrect contours, it is wise to create a list
            # of all possible lines that are related to a new point, then to
            # pick the one that is closest, or, alternatively, create a new
            # seperate line with the intersection as an intial point. For now I
            # will keep the code simple to reduce unneccesary work
            for iNewLower in range(0, len(newLowerCombined[iSection])):
                found = False

                for iLowerLine in range(0, len(linesLower)):
                    if abs(linesLower[iLowerLine][-1][2] - newLowerCombined[iSection][iNewLower][2]) <= 1 and \
                            abs(linesLower[iLowerLine][-1][0] - newLowerCombined[iSection][iNewLower][0]) <= 2:
                        # Point is close enought to other point
                        linesLower[iLowerLine].append(newLowerCombined[iSection][iNewLower])
                        found = True
                        break

                if not found:
                    # Create a new line
                    linesLower.append([newLowerCombined[iSection][iNewLower]])

            for iNewUpper in range(0, len(newUpperCombined[iSection])):
                found = False;

                for iUpperLine in range(0, len(linesUpper)):
                    if abs(linesUpper[iUpperLine][-1][2] - newUpperCombined[iSection][iNewUpper][2]) <= 1 and \
                            abs(linesUpper[iUpperLine][-1][0] - newUpperCombined[iSection][iNewUpper][0]) <= 2:
                        # Point is close enough
                        linesUpper[iUpperLine].append(newUpperCombined[iSection][iNewUpper])
                        found = True
                        break

                if not found:
                    # Create a new upper line
                    linesUpper.append([newUpperCombined[iSection][iNewUpper]])

            # Extend the extents array by repeating the same process as
            # performed after obtaining the initial points
            if len(newLowerCombined[iSection] != 0 or len(newUpperCombined[iSection] != 0)):
                # Create a new extent entry
                newExtent = [height[iSection], []]

                if len(newLowerCombined[iSection]) != len(newUpperCombined[iSection]):
                    if abs(len(newLowerCombined[iSection]) - len(newUpperCombined[iSection])) > 1:
                        raise RuntimeError("Expected length of 'newLowerCombined' to differ " +
                            "by a maximum of one from 'newUpperCombined'")

                    if len(newLowerCombined[iSection]) < len(newUpperCombined[iSection]):
                        # There is one upper-section at the left-most side
                        newExtent[1].append([axisDeltaV[0], newUpperCombined[iSection]])

                        for i in range(0, len(newLowerCombined[iSection])):
                            newExtent[1].append([newLowerCombined[iSection][i][1],
                                                 newUpperCombined[iSection][i + 1][1])
                    else:
                        # There is one lower-section at the right-most side
                        for i in range(0, len(newUpperCombined[iSection])):
                            newExtent[1].append([newLowerCombined[iSection][i][1],
                                                [newUpperCombined[iSection][i][1]]])

                        newExtent[1].append([newLowerCombined[iSection][-1][1], axisDeltaV[-1]])
                else:
                    if newUpperCombined[iSection][0][1] > newLowerCombined[iSection][0][1]:
                        # Default ordering, only sets of lower and
                        # upper-sections exist
                        for i in range(0, len(newUpperCombined[iSection]):
                            newExtent[1].append([newLowerCombined[iSection][i][1],
                                                 newUpperCombined[iSection][i][1]])
                    else:
                        # There is one upper-section to the left-most side and
                        # one upper-section to the right-most side
                        newExtent[1].append([axisDeltaV[0], newUpperCombined[iSection][0][1]])

                        for i in range(1, len(newLowerCombined[iSection])):
                            newExtent[1].append([newLowerCombined[iSection][i - 1][1],
                                                 newUpperCombined[iSection][i][1]])

                        newExtent[1].append([newLowerCombined[iSection][-1][1], axisDeltaV[-1]])

                extents.append(newExtent)

    # Return the results
    return linesLower, linesUpper, extents

def GenerateContour(axisHeight, axisDeltaV, maps, qInfMin, qInfMax,
        alphaMin, alphaMax, PReqMin, PReqMax, vInfOveraMin, vInfOveraMax):
    # Loop through the various maps and find the limits to generate contours for
    # each of the limits. When all contours are generated then combine them into
    # a single contour where all the conditions are satisfied
    return True

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

def __testPlotMaps__(severity):
    atm = Atmosphere.Atmosphere()
    axisHeight = np.linspace(10, 75, 250) * 1e3
    axisDeltaV = np.linspace(-50, 50, 250)
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

def __testDetermineSingleContour__():
    axisHeight = np.linspace(0, 5, 5)
    axisDeltaV = np.linspace(5, 10, 5)
    data = np.asarray([
        [ 1, 5, 3, 5, 1],
        [ 1, 4, 4, 3, 2],
        [ 4, 3, 4, 4, 1],
        [ 2, 3, 2, 1, 1],
        [ 1, 2, 4, 1, 1]
    ])

    linesLower, linesUpper, extents = __determineSingleContour__(axisHeight,
        axisDeltaV, data, 3.5, 6.0)

__testPlotMaps__(-1.5)
__testPlotMaps__(1.5)
__testPlotMaps3D__()
