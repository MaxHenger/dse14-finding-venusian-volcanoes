# -*- coding: utf-8 -*-
"""
Created on Tue May 31 12:06:40 2016

This was one of the initial files used to investigate the valid range of flight
paths the aircraft could take throughout the atmosphere. This file is also the
source of one complete day of worrying, being sad and cussing.

@author: MaxHenger
"""

import Atmosphere

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import TrackLookup
import TrackCommon

def plotMap(axMap, axBar, map, cmapMap, cmapBar, xAxis, yAxis, xLabel, yLabel, zLabel, contours,
			forceMapMin=None, forceMapMax=None):
	# Retrieve minimum and maximum values and construct a normalization
	xMin = min(xAxis)
	xMax = max(xAxis)
	yMin = min(yAxis)
	yMax = max(yAxis)
	zMin = np.min(map)
	zMax = np.max(map)

	if forceMapMin != None:
		zMin = forceMapMin

	if forceMapMax != None:
		zMax = forceMapMax

	norm = mpl.colors.Normalize(forceMapMin, forceMapMax)

	# Show the heatmap
	axMap.imshow(map, extent=[xMin, xMax, yMin, yMax], origin='lower',
		aspect='auto', cmap=cmapMap, norm=norm)
	axMap.set_xlabel(xLabel)
	axMap.set_ylabel(yLabel)
	axMap.grid(True)

	# Determine the colorbar's colors
	colors = []

	if isinstance(cmapBar, str):
		colors = ['k'] * len(contours)
	else:
		for i in contours:
			colors.append(cmapBar((i - zMin) / (zMax - zMin)))

	hContour = axMap.contour(xAxis, yAxis, map, contours, colors=colors)
	axMap.clabel(hContour)

	# Add the colormap
	if axBar != None:
		cbb = mpl.colorbar.ColorbarBase(axBar, cmap=cmapMap, norm=norm)
		cbb.set_label(zLabel, rotation=90, fontsize=15)
		cbb.add_lines(contours, colors, 1.0)

def InvestigateCruise(W, S, lookupCl, lookupCd, inclination):
	atm = Atmosphere.Atmosphere()
	alphaMin = lookupCl.getPoints()[0][0]
	alphaMax = lookupCl.getPoints()[0][-1]

	# Set the settings for graphing
	numPointsHeight = 75
	numPointsDeltaV = 50
	height = np.linspace(25000, 75000, numPointsHeight)
	deltaV = np.linspace(-40.0, 40.0, numPointsDeltaV)

	# Preallocate arrays to store results in
	qInf = np.zeros([numPointsHeight, numPointsDeltaV])
	vInf = np.zeros([numPointsHeight, numPointsDeltaV])
	pReq = np.zeros([numPointsHeight, numPointsDeltaV])
	alpha = np.zeros([numPointsHeight, numPointsDeltaV])

	# Precalculate the velocity and the density at a given height
	vZonal = atm.velocityZonal(height, 0, 0)[1]
	vSpeedOfSound = atm.speedOfSound(height, 0, 0)
	print('speed of sound', vSpeedOfSound)
	rho = atm.density(height, 0, 0)[1]

	for iHeight in range(0, numPointsHeight):
		for iDeltaV in range(0, numPointsDeltaV):
			vCompound = deltaV[iDeltaV] + vZonal[iHeight]
			curQInf = 0.5 * rho[iHeight] * vCompound**2.0
			qInf[iHeight, iDeltaV] = curQInf
			vInf[iHeight, iDeltaV] = vCompound / vSpeedOfSound[iHeight]

			itAlphaMin = alphaMin
			itAlphaMax = alphaMax
			itAlphaRange = itAlphaMax - itAlphaMin
			itAlpha = 0
			alphaBest = 0.0

			for i in range(0, 5):
				alphaTest = np.linspace(itAlphaMin, itAlphaMax, 10)
				alphaDeltaBest = 1e9

				for j in range(0, len(alphaTest)):
					itAlphaNew = (np.arctan2(W / (curQInf * S) - lookupCl(alphaTest[j]),
						lookupCd(alphaTest[j])) - inclination) * 180.0 / np.pi

					if abs(itAlphaNew - alphaTest[j]) < alphaDeltaBest:
						alphaBest = alphaTest[j]
						alphaDeltaBest = abs(itAlphaNew - alphaTest[j])

				itAlphaMin = max(alphaMin, alphaBest - 0.25 * itAlphaRange)
				itAlphaMax = min(alphaMax, alphaBest + 0.25 * itAlphaRange)
				itAlphaRange /= 4

			alpha[iHeight, iDeltaV] = alphaBest
			pReq[iHeight, iDeltaV] = vCompound * curQInf * S * lookupCd(alphaBest) / np.cos(alphaBest / 180.0 * np.pi + inclination)

	# settings for plotting
	border = 0.08
	xdist = 0.08
	ydist = 0.05
	cbardist = 0.02
	cbarwidth = 0.015

	# plot
	plotWidth = 0.5 - border - xdist / 2.0 - cbarwidth - cbardist
	plotHeight = 0.5 - border - ydist / 2.0
	x2 = border + plotWidth + cbarwidth + cbardist + xdist
	y2 = border + plotHeight + ydist

	fig = plt.figure()
	axQInf = fig.add_axes([border, y2, plotWidth, plotHeight])
	axQInfMap = fig.add_axes([border + plotWidth + cbardist / 2, y2, cbarwidth, plotHeight])
	plotMap(axQInf, axQInfMap, np.log10(qInf), plt.get_cmap('jet'),
	 	'k', deltaV, height / 1e3, 'deltaV [m/s]',
		'height [km]', 'qInf [Pa]', [2, 2.5, 3.0, 3.5, 4.0, 4.5], forceMapMin=0.0)

	axPReq = fig.add_axes([x2, y2, plotWidth, plotHeight])
	axPReqMap = fig.add_axes([x2 + plotWidth + cbardist / 2, y2, cbarwidth, plotHeight])
	plotMap(axPReq, axPReqMap, pReq / 1e3, plt.get_cmap('jet'), 'k', deltaV,
		height / 1e3, 'deltaV [m/s]', 'height [km]', 'PReq [kW]', [5, 10, 15, 20, 30, 40, 50, 60],
            forceMapMin=0.0, forceMapMax=100.0)

	axAlpha = fig.add_axes([border, border, plotWidth, plotHeight])
	axAlphaMap = fig.add_axes([border + plotWidth + cbardist / 2, border, cbarwidth, plotHeight])
	plotMap(axAlpha, axAlphaMap, alpha, plt.get_cmap('jet'), 'k', deltaV,
		height / 1e3, 'deltaV [deg]', 'height [km]', 'alpha [deg]',
		[-7.5, -5.0, -2.5, 0, 2.5, 5.0, 7.5],
		forceMapMin=-1.0, forceMapMax=10.0)

	axVInf = fig.add_axes([x2, border, plotWidth, plotHeight])
	axVInfMap = fig.add_axes([x2 + plotWidth + cbardist / 2.0, border, cbarwidth, plotHeight])
	plotMap(axVInf, axVInfMap, vInf, plt.get_cmap('jet'), 'k', deltaV,
		height / 1e3, 'deltaV [deg]', 'height [km]', 'vInf [m/s]',
		[0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])

InvestigateCruise(700*8.8, 35, *TrackCommon.LoadAerodynamicData('./data/aerodynamicPerformance/Cl.csv',
				  './data/aerodynamicPerformance/Cd.csv'), 0)
