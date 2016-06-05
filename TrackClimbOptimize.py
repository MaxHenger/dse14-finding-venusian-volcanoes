# -*- coding: utf-8 -*-
"""
Created on Mon May 30 02:54:42 2016

This file contains the climbing optimization routine for the aircraft. With an
intended lower height, upper height and velocities this routine will attempt to
maximize rate of climb using the given aerodynamic performance while staying
inside provided horizontal velocity bounds.

The intention of this file is to experiment with various different power values
between different ranges of altitudes to figure out to most efficient way to
climb to the desired altitude while minimizing power use.

The results can be stored to a file once the algorithm is done running. It will
write the results using an instance of the 'TrackStorage' class.

@author: MaxHenger
"""

import matplotlib.pyplot as plt
import numpy as np

import Atmosphere
import TrackClimb
import TrackCommon
import TrackBiasMap
import TrackLookup
import TrackStorage
import TrackAngleOfAttack

def OptimizeClimb(heightLower, heightUpper, heightQuit, vHorInitial, vVerInitial,
				  longitude, latitude, W, S, PRequired, inclination, dt, lookupCl,
				  lookupCd, lookupBoundLowerVInf=None, lookupBoundUpperVInf=None,
				  storeResults=True):
	# Construct lookup tables and interpolators
	atmosphere = Atmosphere.Atmosphere()
	lookupdCldAlpha = lookupCl.getDerivative()
	lookupdCddAlpha = lookupCd.getDerivative()
	lookupReverseCl = TrackLookup.LookupSegmented1D(lookupCl.getPoints()[1],
		lookupCl.getPoints()[0])

	# Retrieve ranges of angles of attack from the lookup tables
	numAlpha = 125
	alphaLimits = lookupCl.getPoints()
	alphaLimits = [alphaLimits[0][0], alphaLimits[0][-1]]

	# Test: find the angle of attack with the maximum L/D
	#ClPoints = np.asarray(lookupCl.getPoints())
	#CdPoints = np.asarray(lookupCd.getPoints())
	#alphaMaxClCd = ClPoints[0][np.argmax(ClPoints[1] / CdPoints[1])]

	# Settings for the optimization routine
	biasLimit = 0.10 # percent
	biasStep = 0.15 # percent
	biasWidth = 5000 # km
	percentSpeedOfSound = 0.75

	alphaDotLimit = 0.5 # degree/s
	gammaDotLimit = 1.0 / 180.0 * np.pi # rad/s
	gammaLimit = np.pi / 2.0 # rad

	updateCount = 35 # number of iterations before printing an update statement
	averageTime = 2.5 # number of seconds to average from the results for the resimulation

	# Set initial values
	initialRho = atmosphere.density(heightLower, latitude, longitude)[1]
	initialZonal = atmosphere.velocityZonal(heightLower, latitude, longitude)[1]
	initialVInf = np.sqrt(np.power(initialZonal + vHorInitial, 2.0) + np.power(vVerInitial, 2.0))
	initialGamma = np.arctan2(vVerInitial, initialZonal + vHorInitial)

	initialAlpha = 0.0

	if abs(initialGamma) < 0.1 / 180.0 * np.pi:
		initialAlpha = TrackAngleOfAttack.AngleOfAttackSteady(W, S,
			0.5 * initialRho * initialVInf**2.0, lookupReverseCl)
	else:
		initialAlpha = TrackAngleOfAttack.AngleOfAttackPowered(W, S,
			0.5 * initialRho * initialVInf**2.0, PRequired / initialVInf,
			initialGamma, inclination, lookupCl, lookupdCldAlpha)

	if initialAlpha[1] == False:
		raise ValueError("Initial angle of attack is invalid")

	initialAlpha = initialAlpha[0]

	alpha = [initialAlpha]
	vHor = [vHorInitial]
	vVer = [vVerInitial]
	height = [heightLower]
	gamma = [initialGamma]
	time = [0.0]
	vLim = [0.0]

	# Set bias maps and their associated variables
	biasBaseDefaultGamma = 0.975
	biasBaseDefaultVInf = 0.975
	biasBaseDefaultGammaDot = 0.975
	biasBaseDefaultBoundVInf = 0.975

	biasBaseGamma = biasBaseDefaultGamma
	biasBaseVInf = biasBaseDefaultVInf
	biasBaseGammaDot = biasBaseDefaultGammaDot
	biasBaseBoundLowerVInf = biasBaseDefaultBoundVInf
	biasBaseBoundUpperVInf = biasBaseDefaultBoundVInf

	heightBiasLower = heightQuit
	heightBiasUpper = heightUpper + (heightLower - heightQuit)

	biasGamma = TrackBiasMap.BiasMap("gamma", heightBiasLower, heightBiasUpper, 1024, biasBaseGamma)
	biasVInf = TrackBiasMap.BiasMap("vInf", heightBiasLower, heightBiasUpper, 1024, biasBaseVInf)
	biasGammaDot = TrackBiasMap.BiasMap("gammaDot", heightBiasLower, heightBiasUpper, 1024, biasBaseGammaDot)
	biasBoundLowerVInf = TrackBiasMap.BiasMap("vInfLower", heightBiasLower, heightBiasUpper, 1024, biasBaseBoundLowerVInf)
	biasBoundUpperVInf = TrackBiasMap.BiasMap("vInfUpper", heightBiasLower, heightBiasUpper, 1024, biasBaseBoundUpperVInf)

	# Start iterating
	solved = False

	while not solved:
		totalTime = 0.0

		alpha = [initialAlpha]
		vHor = [vHorInitial]
		vVer = [vVerInitial]
		height = [heightLower]
		gamma = [initialGamma]
		time = [0.0]
		vLim = [0.0]

		gammaOld = gamma[0]

		iIteration = 0
		solved = True

		while height[-1] < heightUpper:
			# Determine the new valid range of angles of attack
			alphaNew = np.linspace(np.max([alphaLimits[0], alpha[-1] - dt * alphaDotLimit]),
								   np.min([alphaLimits[1], alpha[-1] + dt * alphaDotLimit]), numAlpha)

			# Determine the new flight variables
			vHorNew, vVerNew, gammaNew, hNew = TrackClimb.Step(height[-1],
				alpha[-1], gamma[-1], vHor[-1], vVer[-1], longitude, latitude,
				PRequired, W, S, inclination, alphaNew, dt, lookupCl, lookupCd,
				atmosphere, tol=1e-8, relax=0.8)

			totalTime = dt * (iIteration + 1)

			# Filter the valid solutions from the invalid ones. Keep track of
			# the number of offenders to know how to alter the bias maps
			iValid = []
			vZonal = atmosphere.velocityZonal(hNew, latitude, longitude)[1]
			vInf = np.sqrt(np.power(vHorNew + vZonal, 2.0) + np.power(vVerNew, 2.0))
			vLimit = atmosphere.speedOfSound(hNew, latitude, longitude) * percentSpeedOfSound

			gammaOffenders = 0
			vInfOffenders = 0
			gammaDotOffenders = 0
			boundLowerVInfOffenders = 0
			boundUpperVInfOffenders = 0

			for i in range(0, len(alphaNew)):
				# Determine if this solution is valid. If any of the conditions
				# fails then take note of it
				isOffender = False

				if gammaNew[i] < -gammaLimit or gammaNew[i] > gammaLimit:
					gammaOffenders += 1
					isOffender = True

				if vInf[i] > vLimit[i]:
					vInfOffenders += 1
					isOffender = True

				if iIteration >= int(5 / dt):
					if abs((gammaNew[i] - gammaOld) / dt) > gammaDotLimit:
						gammaDotOffenders += 1
						isOffender = True

				if lookupBoundLowerVInf != None:
					if vHorNew[i] < lookupBoundLowerVInf(hNew[i]):
						boundLowerVInfOffenders += 1
						isOffender = True

				if lookupBoundUpperVInf != None:
					if vHorNew[i] > lookupBoundUpperVInf(hNew[i]):
						boundUpperVInfOffenders += 1
						isOffender = True

				if isOffender:
					continue

				# This item is valid
				iValid.append(i)

			if len(iValid) == 0:
				# No valid solutions were found
				print("\n * Did not find a solution at:\n > " +
					TrackCommon.StringPad("t = ", totalTime, 3, 10) + " s \n > " +
					TrackCommon.StringPad("h = ", hNew[-1], 3, 10) + ' m\n')

				toShow = int(len(alphaNew) / 2)
				print(TrackCommon.StringPad(" > gamma          = ", gammaNew[toShow] * 180.0 / np.pi, 3, 10) + " deg")
				print(TrackCommon.StringPad(" > vInf           = ", vInf[toShow], 3, 10) + " m/s")
				print(TrackCommon.StringPad(" > vLimit         = ", vLimit[toShow], 3, 10) + " m/s")
				print(TrackCommon.StringPad(" > vZonal         = ", vZonal[toShow], 3, 10) + " m/s")
				print(TrackCommon.StringPad(" > vHor           = ", vHorNew[toShow], 3, 10) + " m/s")
				print(TrackCommon.StringPad(" > vVer           = ", vVerNew[toShow], 3, 10) + " m/s")
				print(TrackCommon.StringPad(" > gammaDot       = ", abs((gammaNew[toShow] - gammaOld) * 180.0 / np.pi / dt), 5, 10) + " deg/s")
				print(TrackCommon.StringPad(" > gammaDot limit = ", gammaDotLimit * 180.0 / np.pi, 5, 10) + " deg/s")

				if lookupBoundLowerVInf != None:
					print(TrackCommon.StringPad(" > vHor lower     = ", lookupBoundLowerVInf(hNew[toShow]), 3, 10) + " m/s")

				if lookupBoundUpperVInf != None:
					print(TrackCommon.StringPad(" > vHor upper     = ", lookupBoundUpperVInf(hNew[toShow]), 3, 10) + " m/s")

				# Determine how to adjust the bias maps
				listOffenders = [gammaOffenders, vInfOffenders, gammaDotOffenders,
					boundLowerVInfOffenders, boundUpperVInfOffenders]
				iWorstOffender = np.argmax(listOffenders)

				if listOffenders[iWorstOffender] == 0:
					print('\n * No offenders, adjusting base biases:')

					biasBaseGamma, biasBaseVInf, biasBaseGammaDot = \
						TrackCommon.AdjustBiasMapCommonly([biasGamma, biasVInf,
						biasGammaDot, biasBoundLowerVInf, biasBoundUpperVInf],
						biasStep, ['gamma', 'vInf', 'gammaDot', 'vLower', 'vUpper'])

					print('')
				else:
					# For the reason behind the following indices, see the
					# listOffenders variable above
					print('\n * Adjusting bias, average gamma:', np.average(abs(gammaNew)) * 180.0 / np.pi)

					if iWorstOffender == 0:
						# Adjust gamma bias locally
						biasBaseGamma = TrackCommon.AdjustBiasMapIndividually(biasGamma,
							biasStep, height[-1], biasWidth, 'gamma')
					elif iWorstOffender == 1:
						# Adjust vInf bias locally
						biasBaseVInf = TrackCommon.AdjustBiasMapIndividually(biasVInf,
							biasStep, height[-1], biasWidth, 'vInf')
					elif iWorstOffender == 2:
						# Adjust gammaDot bias locally
						biasBaseGammaDot = TrackCommon.AdjustBiasMapIndividually(biasGammaDot,
							biasStep, height[-1], biasWidth, 'gammaDot')
					elif iWorstOffender == 3:
						# Adjust lower vInf bias
						biasBaseBoundLowerVInf = TrackCommon.AdjustBiasMapIndividually(
							biasBoundLowerVInf, biasStep, height[-1], biasWidth, 'vLower')
					elif iWorstOffender == 4:
						# Adjust upper vInf bias
						biasBaseBoundUpperVInf = TrackCommon.AdjustBiasMapIndividually(
							biasBoundUpperVInf, biasStep, height[-1], biasWidth, 'vUpper')
					else:
						raise RuntimeError("Unrecognized offender index for bias map")

					print('')

				if biasBaseGamma < biasLimit or biasBaseVInf < biasLimit or \
					biasBaseGammaDot < biasLimit or biasBaseBoundLowerVInf < biasLimit or \
					biasBaseBoundUpperVInf < biasLimit:
					print("Failed to find a solution")
					break

				# Restart with a new bias
				solved = False
				break

			# In the valid solutions maximize a certain metric (very subject to
			# change. This is screwing around, not math!)
			bestMetric = 1e19
			iSolution = 0

			for i in iValid:
				# New attempt at constructing a metric:
				# - base metric: maximizing positive vertical speed (which, for
				# 	a constant deltaTime, is the same as maximizing deltaHeight)
				metric = ((vLimit[i] - vVerNew[i]) / vLimit[i])**2.0
				#metric = 0.0

				# - influence of dgamma/dt
				dgammadt = (gammaNew[i] - gammaOld) / dt
				curBiasGammaDot = biasGammaDot(hNew[i])
				if dgammadt > curBiasGammaDot * gammaDotLimit:
					metric += ((dgammadt - gammaDotLimit * curBiasGammaDot) /
						(gammaDotLimit * (1.0 - curBiasGammaDot)))**2.0

				# - influence of freestream velocity's proximity to the limit
				curBiasVInf = biasVInf(hNew[i])
				if vInf[i] > vLimit[i] * curBiasVInf:
					metric += ((vInf[i] - vLimit[i] * curBiasVInf) /
						(vLimit[i] * (1.0 - curBiasGammaDot)))**2.0

				# - influence of flight path angle's proximity to the limit
				curBiasGamma = biasGamma(hNew[i])
				if abs(gammaNew[i]) > curBiasGamma * gammaLimit:
					metric += ((abs(gammaNew[i]) - curBiasGamma * gammaLimit) /
						(gammaLimit * (1 - curBiasGamma)))**2.0

				# - influence of the lower vInf bound
				if lookupBoundLowerVInf != None:
					curBiasBoundLowerVInf = biasBoundLowerVInf(hNew[i])
					curBoundLowerVInf = lookupBoundLowerVInf(hNew[i])

					divFactor = vLimit[i] - curBoundLowerVInf

					if lookupBoundUpperVInf != None:
						divFactor = lookupBoundUpperVInf(hNew[i]) - curBoundLowerVInf

					curBoundDelta = (1.0 - curBiasBoundLowerVInf) * divFactor

					if vHorNew[i] - curBoundLowerVInf < curBoundDelta:
						metric += ((curBoundLowerVInf + curBoundDelta - vHorNew[i]) / curBoundDelta)**2.0

				# - influence of the upper vInf bound
				if lookupBoundUpperVInf != None:
					curBiasBoundUpperVInf = biasBoundUpperVInf(hNew[i])
					curBoundUpperVInf = lookupBoundUpperVInf(hNew[i])

					divFactor = curBoundUpperVInf

					if lookupBoundLowerVInf != None:
						divFactor = curBoundUpperVInf - lookupBoundLowerVInf(hNew[i])

					curBoundDelta = (1.0 - curBiasBoundUpperVInf) * divFactor

					if curBoundUpperVInf - vHorNew[i] < curBoundDelta:
						#print('upper: ' +
						#	TrackCommon.StringPad(' hor: ', vHorNew[i], 1, 6) + " m/s" +
						#	TrackCommon.StringPad(', lim: ', curBoundUpperVInf, 1, 6) + " m/s" +
						#	TrackCommon.StringPad(', div: ', divFactor, 1, 6) + " m/s" +
						#	TrackCommon.StringPad(', metric: ', ((vHorNew[i] - curBoundUpperVInf) / divFactor - curBiasBoundUpperVInf)**2.0, 3, 6))

						metric += ((vHorNew[i] + curBoundDelta - curBoundUpperVInf) / curBoundDelta)**2.0

				# TODO: TEST DISTANCE FROM ALPHA CL/CD MAX
				#metric += ((alphaNew[i] - alphaMaxClCd) / ClPoints[0][-1] *
				#		   vInf[i] / vLimit[i])**2.0

				if metric < bestMetric:
					iSolution = i
					bestMetric = metric

			# Append the new values to the solution arrays
			height.append(hNew[iSolution])
			alpha.append(alphaNew[iSolution])
			gamma.append(gammaNew[iSolution])
			vHor.append(vHorNew[iSolution])
			vVer.append(vVerNew[iSolution])
			vLim.append(vLimit[iSolution])
			time.append(totalTime)

			gammaOld = gammaNew[iSolution]

			if iIteration % updateCount == 0:
				print(TrackCommon.StringPad("Solved at t = ", totalTime, 3, 8) +
					  TrackCommon.StringPad(" s, h = ", hNew[iSolution], 0, 7) +
					  TrackCommon.StringPad(" m, Vver = ", vVerNew[iSolution], 2, 6) +
					  TrackCommon.StringPad(" m/s, gamma = ", gammaNew[iSolution] * 180.0 / np.pi, 3, 8) +
					  TrackCommon.StringPad(" deg, alpha = ", alphaNew[iSolution], 3, 8) + " deg")

			iIteration += 1

			# If this position is reached then the algorithm iterated through
			# all the height values. The 'solved' variable will take care of
			# quitting the main processing loop or not

	# Rerun the simulation with values that are averaged over a prolonged period
	# of time to see if the results match
	alphaFinal = [alpha[0]]
	vHorFinal = [vHor[0]]
	vVerFinal = [vVer[0]]
	heightFinal = [height[0]]
	gammaFinal = [gamma[0]]
	timeFinal = [time[0]]
	vLimFinal = [vLim[0]]

	if True:
		numSteps = int(averageTime / dt)

		for iIteration in range(0, len(alpha)):
			# Average out angles of attack
			alphaAverage = 0.0
			alphaMin = int(max(iIteration - numSteps, 0))
			alphaMax = int(min(iIteration + numSteps, len(alpha)))

			for i in range(alphaMin, alphaMax):
				alphaAverage += alpha[i]

			alphaAverage /= (alphaMax - alphaMin)

			# Perform a step in the simulation
			vHorNew, vVerNew, gammaNew, hNew = TrackClimb.Step(heightFinal[-1],
				alphaFinal[-1], gammaFinal[-1], vHorFinal[-1], vVerFinal[-1], longitude,
				latitude, PRequired, W, S, inclination, np.asarray([alphaFinal[-1]]),
				dt, lookupCl, lookupCd, atmosphere, tol=1e-8, relax=0.8)

			# Store results
			alphaFinal.append(alphaAverage)
			vHorFinal.append(vHorNew[0])
			vVerFinal.append(vVerNew[0])
			heightFinal.append(hNew[0])
			gammaFinal.append(gammaNew[0])
			timeFinal.append(time[iIteration])

	# Plot the results
	fig = plt.figure()
	axAlpha = fig.add_subplot(321)
	axGamma = fig.add_subplot(322)
	axHeight = fig.add_subplot(323)
	axVer = fig.add_subplot(324)
	axHor = fig.add_subplot(325)
	axVinf = fig.add_subplot(326)

	axAlpha.plot(time, alpha, 'g')
	axAlpha.plot(timeFinal, alphaFinal, 'r')
	axAlpha.set_xlabel('time [s]')
	axAlpha.set_ylabel('alpha [deg]')
	axAlpha.grid(True)

	axGamma.plot(time, np.asarray(gamma) * 180.0 / np.pi, 'g')
	axGamma.plot(timeFinal, np.asarray(gammaFinal) * 180.0 / np.pi, 'r')
	axGamma.set_xlabel('time [s]')
	axGamma.set_ylabel('gamma [deg]')
	axGamma.grid(True)

	axHeight.plot(time, np.asarray(height) / 1e3, 'g')
	axHeight.plot(timeFinal, np.asarray(heightFinal) / 1e3, 'r')
	axHeight.set_xlabel('time [s]')
	axHeight.set_ylabel('height [km]')
	axHeight.grid(True)

	axVer.plot(time, vVer, 'g')
	axVer.plot(timeFinal, vVerFinal, 'r')
	axVer.set_xlabel('time [s]')
	axVer.set_ylabel('vertical speed [m/s]')
	axVer.grid(True)

	# prepare vHor
	boundLowerVInf = np.zeros([len(height)])
	boundUpperVInf = np.zeros([len(height)])

	for i in range(0, len(height)):
		boundLowerVInf[i] = lookupBoundLowerVInf(height[i])
		boundUpperVInf[i] = lookupBoundUpperVInf(height[i])

	lVHor, = axHor.plot(time, vHor, 'g')
	lVHorFinal, = axHor.plot(timeFinal, vHorFinal, 'r')
	lVLower, = axHor.plot(time, boundLowerVInf, 'k--')
	lVUpper, = axHor.plot(time, boundUpperVInf, 'k--')
	axHor.set_xlabel('time [s]')
	axHor.set_ylabel('horizontal speed [m/s]')
	axHor.grid(True)

	# prepare Vinf
	vInf = np.zeros([len(vVer)])

	for i in range(0, len(vVer)):
		vZonal = atmosphere.velocityZonal(height[i], latitude, longitude)[1]
		vInf[i] = np.sqrt(np.power(vZonal + vHor[i], 2.0) + np.power(vVer[i], 2.0))

	vInfFinal = np.zeros([len(vVerFinal)])

	for i in range(0, len(vVerFinal)):
		vZonal = atmosphere.velocityZonal(heightFinal[i], latitude, longitude)[1]
		vInfFinal[i] = np.sqrt(np.power(vZonal + vHorFinal[i], 2.0) + np.power(vVerFinal[i], 2.0))

	axVinf.plot(time, vInf, 'g')
	axVinf.plot(timeFinal, vInfFinal, 'r')
	axVinf.plot(time, vLim, 'g--')
	axVinf.grid(True)

	fig.suptitle("Preq = " + str(round(PRequired / 1e3, 3)) + " kW")

	# Plot the bias maps
	fig = plt.figure()
	axBias = fig.add_subplot(111)
	lGamma, = axBias.plot(biasGamma.getAxis() / 1e3, biasGamma.getMap() * 1e2, 'r')
	lVInf, = axBias.plot(biasVInf.getAxis() / 1e3, biasVInf.getMap() * 1e2, 'g')
	lGammaDot, = axBias.plot(biasGammaDot.getAxis() / 1e3, biasGammaDot.getMap() * 1e2, 'b')
	lBoundLowerVInf, = axBias.plot(biasBoundLowerVInf.getAxis() / 1e3, biasBoundLowerVInf.getMap() * 1e2, 'y')
	lBoundUpperVInf, = axBias.plot(biasBoundUpperVInf.getAxis() / 1e3, biasBoundUpperVInf.getMap() * 1e2, 'k')

	axBias.legend([lGamma, lVInf, lGammaDot, lBoundLowerVInf, lBoundUpperVInf],
				  ['gamma', 'vInf', 'gammaDot', 'vLower', 'vUpper'])
	axBias.set_xlabel('h [km]')
	axBias.set_ylabel('bias [%]')
	axBias.grid(True)

	if storeResults:
		file = TrackStorage.DataStorage()
		file.addVariable('time', time)
		file.addVariable('alpha', alpha)
		file.addVariable('gamma', gamma)
		file.addVariable('height', height)
		file.addVariable('vVer', vVer)
		file.addVariable('vHor', vHor)
		file.addVariable('vInf', vInf)
		file.addVariable('biasGamma', biasGamma.getMap(), [biasGamma.getAxis()])
		file.addVariable('biasVInf', biasVInf.getMap(), [biasVInf.getAxis()])
		file.addVariable('biasGammaDot', biasGammaDot.getMap(), [biasGammaDot.getMap()])
		file.addVariable('biasBoundLowerVInf', biasBoundLowerVInf.getMap(), [biasBoundLowerVInf.getAxis()])
		file.addVariable('biasBoundUpperVInf', biasBoundUpperVInf.getMap(), [biasBoundUpperVInf.getAxis()])
		file.addVariable('dt', dt)
		file.addVariable('vLim', vLim)
		file.addVariable('updateCount', updateCount)
		file.addVariable('averageTime', averageTime)
		file.save('climb_' + str(heightLower) + 'to' + str(heightUpper) +
			'_' + str(PRequired) + '_' + str(vHorInitial) +
			'_' + str(vVerInitial) + '.dat')

def PlotClimb(filename):
	file = TrackStorage.DataStorage()
	file.load(filename)
	time = file.getVariable("time").getValues()
	alpha = file.getVariable("alpha").getValues()
	gamma = file.getVariable("gamma").getValues()
	height = file.getVariable("height").getValues()
	vVer = file.getVariable("vVer").getValues()
	vHor = file.getVariable("vHor").getValues()
	vInf = file.getVariable("vInf").getValues()
	vLimit = file.getVariable("vLim").getValues()
	biasAxis = file.getVariable("biasGamma").getAxis(0)
	biasGamma = file.getVariable("biasGamma").getValues()
	biasVInf = file.getVariable("biasVInf").getValues()
	biasGammaDot = file.getVariable("biasGammaDot").getValues()
	biasBoundLowerVInf = file.getVariable("biasBoundLowerVInf").getValues()
	biasBoundUpperVInf = file.getVariable("biasBoundUpperVInf").getValues()
	dt = file.getVariable('dt').getValues()

	fig = plt.figure()
	axSpeedDir = fig.add_subplot(221)
	axSpeed = axSpeedDir.twinx()
	axHeight = fig.add_subplot(222)
	axAnglesLeft = fig.add_subplot(223)
	axAnglesRight = axAnglesLeft.twinx()
	axAnglesDot = fig.add_subplot(224)

	# Plot all relevant speeds
	lSpeedVer, = axSpeed.plot(time, vVer, 'r', label=r'$V_{\mathrm{ver}}$')
	lSpeedHor, = axSpeed.plot(time, vHor, 'g', label=r'$V_{\mathrm{hor}}$')
	lSpeedInf, = axSpeed.plot(time, vInf, 'b', label=r'$V_{\mathrm{\infty}}$')
	lSpeedLim, = axSpeed.plot(time, vLimit, 'k--', label=r'$V_{\mathrm{lim}}$')
	axSpeed.set_xlabel(r'$t\;[s]$')
	axSpeed.set_ylabel(r'$V\;[m/s]$')
	axSpeed.grid(True)
	axSpeed.legend()

	# Plot the height
	axHeight.plot(time, height / 1e3, 'r')
	axHeight.set_xlabel(r'$t\;[s]$')
	axHeight.set_ylabel(r'$h\;[km]$')
	axHeight.grid(True)

	# Plot the angles
	axAnglesLeft.plot(time, alpha, 'r', label=r'$\alpha$')
	for tick in axAnglesLeft.get_yticklabels():
		tick.set_color('r')

	axAnglesRight.plot(time, gamma, 'g', label=r'$\gamma$')
	for tick in axAnglesRight.get_yticklabels():
		tick.set_color('g')

	axAnglesLeft.set_xlabel(r'$t\;[s]$')
	axAnglesLeft.set_ylabel(r'$\alpha\;[\degree]$')
	axAnglesRight.set_ylabel(r'$\gamma\;[\degree]$')

	# Plot the delta angles
	alphaDot = (alpha[1:] - alpha[0:-1]) / dt
	gammaDot = (gamma[1:] - gamma[0:-1]) / dt * 180.0 / np.pi
	axAnglesDot.plot(time[0:-1], alphaDot, 'r', label=r'$\mathrm{d}\alpha/\mathrm{d}t$')
	axAnglesDot.plot(time[0:-1], gammaDot, 'g', label=r'$\mathrm{d}\gamma/\mathrm{d}t$')
	axAnglesDot.set_xlabel(r'$t\;[s]$')
	axAnglesDot.set_ylabel(r'$\omega\;[\degree/s]$')
	axAnglesDot.legend()

	# Plot the bias maps
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(biasAxis, biasGamma, 'r', label=r'$b_\mathrm{\gamma}$')
	ax.plot(biasAxis, biasVInf, 'g', label=r'$b_\mathrm{V_\infty}$')
	ax.plot(biasAxis, biasGammaDot, 'b', label=r'$b_\mathrm{\dot{\gamma}}$')
	ax.plot(biasAxis, biasBoundLowerVInf, 'y', label=r'$b_{V_{\infty,\mathrm{low}}}$')
	ax.plot(biasAxis, biasBoundUpperVInf, 'k', label=r'$b_{V_{\infty,\mathrm{high}}}$')
	ax.set_xlabel(r'$h\;[km]$')
	ax.set_ylabel(r'$b\;[\%]$')

def __TestOptimizeClimb__():
	lookupCl, lookupCd = TrackCommon.LoadAerodynamicData("./data/aerodynamicPerformance/Cl.csv",
														 "./data/aerodynamicPerformance/Cd.csv")
	lookupLower = TrackLookup.Lookup1D([30000, 47000, 64000], [-15.0, -20.0, -5.0])
	lookupUpper = TrackLookup.Lookup1D([30000, 40000, 50000, 64000], [-5.0, -2.5, -4.0, 30.0])
	#OptimizeClimb(38000, 44000, 30000, -10, 0, 0, 0, 700*8.8, 35.0, 10000, 0,
	#	0.25, lookupCl, lookupCd, lookupLower, lookupUpper)
	OptimizeClimb(38000, 62000, 30000, -10, 0, 0, 0, 700*8.8, 35.0, 40000, 0,
		0.25, lookupCl, lookupCd, lookupLower, lookupUpper)
	'''OptimizeClimb(38000, 44000, 30000, -10, 0, 0, 0, 700*8.8, 35.0, 30000, 0,
		0.25, lookupCl, lookupCd, lookupLower, lookupUpper)
	OptimizeClimb(38000, 44000, 30000, -10, 0, 0, 0, 700*8.8, 35.0, 40000, 0,
		0.25, lookupCl, lookupCd, lookupLower, lookupUpper)
	OptimizeClimb(38000, 44000, 30000, -10, 0, 0, 0, 700*8.8, 35.0, 50000, 0,
		0.25, lookupCl, lookupCd, lookupLower, lookupUpper)
	OptimizeClimb(38000, 44000, 30000, -10, 0, 0, 0, 700*8.8, 35.0, 60000, 0,
		0.25, lookupCl, lookupCd, lookupLower, lookupUpper)
	OptimizeClimb(38000, 44000, 30000, -10, 0, 0, 0, 700*8.8, 35.0, 70000, 0,
		0.25, lookupCl, lookupCd, lookupLower, lookupUpper)'''

#__TestOptimizeClimb__()
#PlotClimb('climb_35000to50000_50000_-10_0.dat')
