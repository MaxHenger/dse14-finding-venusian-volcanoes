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

import matplotlib.pyplot as plt
import matplotlib as mpl
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
              cmap='jet', contours=None, forceNormMin=None, forceNormMax=None):
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
				 lookupCl, lookupCd, atm):
    # Preallocate the resulting maps
    qInf = np.zeros([len(axisHeight), len(axisDeltaV)])
	vInf = np.zeros(qInf.shape)
	alpha = np.zeros(qInf.shape)
	Preq = np.zeros(qInf.shape)

    velocityZonal = atm.velocityZonal(axisHeight, latitude, longitude)
    density = atm.density(axisHeight, latitude, longitude)

    # Adjust the velocity map with respect to the 'relative severity'
    if relativeSeverity >= 0.0:
        velocityZonal = velocityZonal[1] + (velocityZonal[2] - velocityZonal[1]) * relativeSeverity
        density = density[1] + (density[2] - density[1]) * relativeSeverity
    else:
        velocityZonal = velocityZonal[1] + (velocityZonal[1] - velocityZonal[0]) * relativeSeverity
        density = density[1] + (density[1] - density[0]) * relativeSeverity

    for iHeight in range(0, len(axisHeight)):
		# Calculate freesteram velocity and dynamic pressure
		vInfCur = velocityZonal[iHeight] + axisDeltaV
		vInf[iHeight, :] = vInfCur
        qInf[iHeight, :] = 0.5 * density[iHeight] * np.power(vInfCur, 2.0)

		# Angle of attack and power required calculations are performed per
		# velocity value
		for iDeltaV in range(0, len(axisDeltaV)):
			alpha[iHeight, iDeltaV], _ = TrackCommon.AngleOfAttack()

    return result
