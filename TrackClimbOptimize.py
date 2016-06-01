# -*- coding: utf-8 -*-
"""
Created on Mon May 30 02:54:42 2016

@author: MaxHenger
"""

import matplotlib.pyplot as plt
import numpy as np

import Atmosphere
import TrackClimb
import TrackCommon
import TrackBiasMap

def OptimizeClimb(heightLower, heightUpper, heightQuit, vHorInitial, vVerInitial,
				  longitude, latitude, W, S, PRequired, inclination, dt, lookupCl,
				  lookupCd, checkpointHeight=[], checkpointVHor=[],
				  checkpointOffset = 5, storeResults=True):
	# Check input and modify checkpoint definitions to reduce computational work
	if len(checkpointHeight) != len(checkpointVHor):
		raise ValueError("Expected 'checkpointHeight' and 'checkpointVHor' to be " +
			"of the same length")

	hasCheckpoints = False

	if len(checkpointHeight) != 0:
		checkpointHeight = np.asarray(checkpointHeight)
		checkpointVHor = np.asarray(checkpointVHor)
		checkpointOrderIndices = np.argsort(checkpointHeight)
		checkpointHeight = checkpointHeight[checkpointOrderIndices]
		checkpointVHor = checkpointVHor[checkpointOrderIndices]
		hasCheckpoints = True

	atmosphere = Atmosphere.Atmosphere()

	# Retrieve ranges of angles of attack from the lookup tables
	numAlpha = 125
	alphaLimits = lookupCl.getPoints()
	alphaLimits = [alphaLimits[0][0], alphaLimits[0][-1]]

	# Test: find the angle of attack with the maximum L/D
	ClPoints = np.asarray(lookupCl.getPoints())
	CdPoints = np.asarray(lookupCd.getPoints())
	alphaMaxClCd = ClPoints[0][np.argmax(ClPoints[1] / CdPoints[1])]

	# Settings for the optimization routine
	biasLimit = 0.10 # percent
	biasStep = 0.25 # percent
	biasWidth = 5000 # km
	percentSpeedOfSound = 0.75

	alphaDotLimit = 0.5 # degree/s
	gammaDotLimit = 1.5 / 180.0 * np.pi # rad/s
	gammaLimit = np.pi / 2.0 # rad

	checkpointProximity = 500

	updateCount = 35 # number of iterations before printing an update statement

	# Set initial values
	initialZonal = atmosphere.velocityZonal(heightLower, latitude, longitude)[1]
	initialGamma = np.arctan2(vVerInitial, initialZonal + vHorInitial)
	alpha = [0.0]
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

	biasBaseGamma = biasBaseDefaultGamma
	biasBaseVInf = biasBaseDefaultVInf
	biasBaseGammaDot = biasBaseDefaultGammaDot
	biasBaseCheckpoint = [0.975] * len(checkpointHeight)

	biasGamma = TrackBiasMap.BiasMap("gamma", heightQuit, heightUpper + (heightLower - heightQuit), 1024, biasBaseGamma)
	biasVInf = TrackBiasMap.BiasMap("vInf", heightQuit, heightUpper + (heightLower - heightQuit), 1024, biasBaseVInf)
	biasGammaDot = TrackBiasMap.BiasMap("gammaDot", heightQuit, heightUpper + (heightLower - heightQuit), 1024, biasBaseGammaDot)

	biasCheckpoints = []

	for i in range(0, len(checkpointHeight)):
		biasCheckpoints.append(TrackBiasMap.BiasMap("bias" + str(round(checkpointHeight[i], 1)),
			heightQuit, heightUpper + (heightLower - heightQuit), 1024, biasBaseCheckpoint[i]))

	# Start iterating
	solved = False

	while not solved:
		totalTime = 0.0

		alpha = [0.0]
		vHor = [vHorInitial]
		vVer = [vVerInitial]
		height = [heightLower]
		gamma = [initialGamma]
		time = [0.0]
		vLim = [0.0]

		gammaOld = gamma[0]
		iFirstCheckpoint = 0

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
			checkpointOffenders = 0
			checkpointAdjustIndices = [False] * len(biasCheckpoints)

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

				if hasCheckpoints:
					# Loop through checkpoints to ascertain the desired speed is reached
					for iCheckpoint in range(iFirstCheckpoint, len(checkpointHeight)):
						if hNew[i] < checkpointHeight[iCheckpoint] - checkpointProximity:
							# Checkpoints are sorted and the new height is below
							# the first checkpoint, stop looking for them
							break
						elif hNew[i] < checkpointHeight[iCheckpoint] + checkpointProximity:
							# Current height lies within checkpoint
							if abs(vHorNew[i] - checkpointVHor[iCheckpoint]) > checkpointOffset:
								checkpointOffenders += 1
								checkpointAdjustIndices[iCheckpoint] = True
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

				if hasCheckpoints:
					for iCheckpoint in range(iFirstCheckpoint, len(checkpointHeight)):
						print(TrackCommon.StringPad(" > vHor check     = ", checkpointVHor[iCheckpoint], 3, 10) + " m/s")

				# Determine how to adjust the bias maps
				listOffenders = [gammaOffenders, vInfOffenders,
					gammaDotOffenders, checkpointOffenders]
				iWorstOffender = np.argmax(listOffenders)

				if listOffenders[iWorstOffender] == 0:
					print('\n * No offenders, adjusting base biases:')
					biasBaseGamma, biasBaseVInf, biasBaseGammaDot = \
						TrackCommon.AdjustBiasMapCommonly([biasGamma, biasVInf, biasGammaDot],
														  biasStep, ['gamma', 'vInf', 'gammaDot'])

					for iCheckpoint in range(0, len(biasCheckpoints)):
						print('adjusting checkpoint', iCheckpoint, 'bias from',
							round(biasBaseCheckpoint[iCheckpoint], 3), 'to',
							round(biasBaseCheckpoint[iCheckpoint] - biasStep, 3))

						biasBaseCheckpoint[iCheckpoint] -= biasStep
						biasCheckpoints.reset(biasBaseCheckpoint[iCheckpoint])

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
						# Adjust checkpoint biases
						for iCheckpoint in range(iFirstCheckpoint, len(biasCheckpoints)):
							if checkpointAdjustIndices[iCheckpoint]:
								biasBaseCheckpoint[iCheckpoint] = \
									TrackCommon.AdjustBiasMapIndividually(biasCheckpoints[iCheckpoint],
									biasStep, height[-1], biasWidth, 'checkpoint ' + str(iCheckpoint))
					else:
						raise RuntimeError("Unrecognized offender index for bias map")

					print('')

				if biasBaseGamma < biasLimit or biasBaseVInf < biasLimit or \
					biasBaseGammaDot < biasLimit or np.min(biasBaseCheckpoint) < biasLimit:
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
				#metric = ((vVerNew[i] - vLimit[i]) / vLimit[i])**2.0
				metric = 0.0

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

				# - influence of any valid checkpoints
				for iCheckpoint in range(iFirstCheckpoint, len(biasCheckpoints)):
					curBiasCheckpoint = biasCheckpoints[iCheckpoint](hNew[i])

					if hNew[i] < checkpointHeight[iCheckpoint] - \
							(1 - curBiasCheckpoint) * \
							(checkpointHeight[iCheckpoint] - heightLower):
						# height is below the region of interest
						break
					elif hNew[i] < checkpointHeight[iCheckpoint] + checkpointProximity:
						# height is inside the region of interest
						metric += ((vHorNew[i] - checkpointVHor[iCheckpoint]) / vInf[i])**2.0

				# TODO: TEST DISTANCE FROM ALPHA CL/CD MAX
				metric += ((alphaNew[i] - alphaMaxClCd) / ClPoints[0][-1] *
						   vInf[i] / vLimit[i])**2.0

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

			# Check if the first checkpoint to check can be incremented
			if hasCheckpoints and height[-1] > checkpointHeight[iFirstCheckpoint] + \
					checkpointProximity and iFirstCheckpoint < len(biasCheckpoints) - 1:
				iFirstCheckpoint += 1

			# If this position is reached then the algorithm iterated through
			# all the height values. The 'solved' variable will take care of
			# quitting the main processing loop or not

	# Plot the results
	fig = plt.figure()
	axAlpha = fig.add_subplot(321)
	axGamma = fig.add_subplot(322)
	axHeight = fig.add_subplot(323)
	axVer = fig.add_subplot(324)
	axHor = fig.add_subplot(325)
	axVinf = fig.add_subplot(326)

	axAlpha.plot(time, alpha, 'g')
	axAlpha.set_xlabel('time [s]')
	axAlpha.set_ylabel('alpha [deg]')
	axAlpha.grid(True)

	axGamma.plot(time, np.asarray(gamma) * 180.0 / np.pi, 'g')
	axGamma.set_xlabel('time [s]')
	axGamma.set_ylabel('gamma [deg]')
	axGamma.grid(True)

	axHeight.plot(time, np.asarray(height) / 1e3, 'g')
	axHeight.set_xlabel('time [s]')
	axHeight.set_ylabel('height [km]')
	axHeight.grid(True)

	axVer.plot(time, vVer, 'g')
	axVer.set_xlabel('time [s]')
	axVer.set_ylabel('vertical speed [m/s]')
	axVer.grid(True)

	axHor.plot(time, vHor, 'g')
	axHor.set_xlabel('time [s]')
	axHor.set_ylabel('horizontal speed [m/s]')
	axHor.grid(True)

	# prepare Vinf
	Vinf = np.zeros([len(vVer)])

	for i in range(0, len(vVer)):
		vZonal = atmosphere.velocityZonal(height[i], latitude, longitude)[1]
		Vinf[i] = np.sqrt(np.power(vZonal + vHor[i], 2.0) + np.power(vVer[i], 2.0))

	axVinf.plot(time, Vinf, 'g')
	axVinf.plot(time, vLim, 'g--')
	axVinf.grid(True)

	# Plot the bias maps
	fig = plt.figure()
	axBias = fig.add_subplot(111)
	lGamma, = axBias.plot(biasGamma.getAxis() / 1e3, biasGamma.getMap() * 1e2, 'r')
	lVInf, = axBias.plot(biasVInf.getAxis() / 1e3, biasVInf.getMap() * 1e2, 'g')
	lGammaDot, = axBias.plot(biasGammaDot.getAxis() / 1e3, biasGammaDot.getMap() * 1e2, 'b')

	lBiases = []
	nBiases = []

	cmap = plt.get_cmap('jet')

	for iCheckpoint in range(0, len(biasCheckpoints)):
		lBias, = axBias.plot(biasCheckpoints[iCheckpoint].getAxis() / 1e3,
			biasCheckpoints[iCheckpoint].getMap() * 1e2,
			color=cmap(iCheckpoint / (len(biasCheckpoints) - 1)))
		lBiases.append(lBias)
		nBiases.append('h = ' + str(round(checkpointHeight[iCheckpoint] / 1e3, 1)))

	lineArray = [lGamma, lVInf, lGammaDot]
	lineArray.extend(lBiases)
	nameArray = ['gamma', 'vInf', 'gammaDot']
	nameArray.extend(nBiases)

	axBias.legend(lineArray, nameArray)
	axBias.set_xlabel('h [km]')
	axBias.set_ylabel('bias [%]')
	axBias.grid(True)

	# TODO: ADD CODE TO SAVE results
	# TODO: ADD CODE TO RESIMULATE AND SMOOTHEN
	# TODO: ADD CODE TO LINEARLY INTERPOLATE BETWEEN CHECKPOINTS

def __TestOptimizeClimb__():
	lookupCl, lookupCd = TrackCommon.LoadAerodynamicData("./data/aerodynamicPerformance/Cl.csv",
														 "./data/aerodynamicPerformance/Cd.csv")
	OptimizeClimb(38000, 44000, 30000, 30, 0, 0, 0, 700*8.8, 35.0, 40000, 0,
		0.25, lookupCl, lookupCd, [39500, 40500], [-20, -15])

__TestOptimizeClimb__()
