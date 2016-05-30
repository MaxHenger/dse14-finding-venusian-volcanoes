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

def OptimizeClimb(heightLower, heightUpper, heightQuit, vHorInitial, vVerInitial,
				  longitude, latitude, W, S, vHorTarget, PRequired, inclination,
				  dt, lookupCl, lookupCd):
	# Retrieve ranges of angles of attack from the lookup tables
	numAlpha = 751
	alphaNew = lookupCl.getPoints()
	alphaNew = np.linspace(alphaNew[0][0], alphaNew[0][-1], numAlpha)
	alphaNew = np.linspace(-5.5, 5.5, numAlpha)

	# Settings for the optimization routine
	biasLimit = 1.0
	biasStep = 0.35
	percentSpeedOfSound = 0.75
	weightVinf = 2.0
	weightGamma = 16.0
	alphaDotLimit = 1.5
	gammaDotLimit = 0.25
	gammaLimit = np.pi / 3.0

	# Set initial values
	alpha = [0.0]
	vHor = [vHorInitial]
	vVer = [vVerInitial]
	height = [heightLower]
	gamma = [0.0]
	time = [0.0]
	vLim = [0.0]

	# Set global variables
	atmosphere = Atmosphere.Atmosphere()

	# Start iterating
	solved = False
	bias = 0.15

	while not solved:
		totalTime = 0.0

		alpha = [0.0]
		vHor = [vHorInitial]
		vVer = [vVerInitial]
		height = [heightLower]
		gamma = [0.0]
		time = [0.0]
		vLim = [0.0]

		alphaOld = alpha[0]
		gammaOld = gamma[0]

		iIteration = 0
		solved = True

		while height[-1] < heightUpper:
			# Determine the new flight variables
			vHorNew, vVerNew, gammaNew, hNew = TrackClimb.Step(height[-1],
				alpha[-1], gamma[-1], vHor[-1], vVer[-1], longitude, latitude,
				PRequired, W, S, inclination, alphaNew, dt, lookupCl, lookupCd,
				atmosphere)

			totalTime = dt * (iIteration + 1)

			# Filter the valid solutions from the invalid ones
			iValid = []
			vZonal = atmosphere.velocityZonal(hNew, latitude, longitude)[1]
			vInf = np.sqrt(np.power(vHorNew + vZonal, 2.0) + np.power(vVerNew, 2.0))
			vLimit = atmosphere.speedOfSound(hNew, latitude, longitude) * percentSpeedOfSound
			vVerMin = 0.0

			for i in range(0, len(alphaNew)):
				# Calculate local freestream velocity and the speed of sound
				if gammaNew[i] > -gammaLimit and gammaNew[i] < gammaLimit and \
						vInf[i] <= vLimit[i] and \
						(iIteration == 0 or (
							abs((alphaNew[i] - alphaOld) / dt) < alphaDotLimit and
							abs((gammaNew[i] - gammaOld) / dt) < gammaDotLimit
						)) and hNew[i] > heightQuit:

					if vVerNew[i] < vVerMin:
						vVerMin = vVerNew[i]

					iValid.append(i)

			if len(iValid) == 0:
				# No valid solutions were found
				bias += biasStep

				print("Did not find a solution at t =", totalTime, ", bias =", bias)
				print(TrackCommon.StringPad("gamma  = ", gammaNew[0] * 180.0 / np.pi, 3, 10) + " deg")
				print(TrackCommon.StringPad("vInf   = ", vInf[0], 3, 10) + " m/s")
				print(TrackCommon.StringPad("vLimit = ", vLimit[0], 3, 10) + " m/s")
				print(TrackCommon.StringPad("vHor   = ", vHor[0], 3, 10) + " m/s")
				print(TrackCommon.StringPad("vVer   = ", vVer[0], 3, 10) + " m/s")

				if bias > biasLimit:
					print("Failed to find a solution")
					break

				# Restart with a new bias
				solved = False
				break

			# In the valid solutions maximize a certain metric (very subject to
			# change. This is screwing around, not math!)
			bestMetric = -1e9
			iSolution = 0

			for i in iValid:
				metric = (vVerNew[i] - vVerMin)**2.0
				metric /= (1 + (0.2*(alphaNew[i] - alphaOld) / dt)**2.0)

				metric /= (1 + (0.1 * ((vVerNew[i] - vVer[-1]) / vVerNew[-1]))**2.0)
				metric /= (1 + (0.1 * ((vHorNew[i] - vHor[-1]) / vHorNew[-1]))**2.0)

				if vInf[i] > vLimit[i] * (1.0 - bias):
					metric *= ((vLimit[i] - vInf[i]) / (bias * vLimit[i]))*2.0
					metric /= weightVinf

				if abs(gammaNew[i]) > gammaLimit * (1.0 - bias):
					metric *= ((gammaLimit - abs(gammaNew[i])) / (bias * gammaLimit))**2.0
					metric /= weightGamma

				if metric > bestMetric:
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
			alphaOld = alphaNew[iSolution]

			if iIteration % 20 == 0:
				print(TrackCommon.StringPad("Solved at t = ", totalTime, 3, 8) +
					  TrackCommon.StringPad(" s, h = ", hNew[iSolution], 0, 7) +
					  TrackCommon.StringPad(" m, Vver = ", vVerNew[iSolution], 2, 6) +
					  TrackCommon.StringPad(" m/s, gamma = ", gammaNew[iSolution] * 180.0 / np.pi, 3, 8) +
					  TrackCommon.StringPad(" deg, alpha = ", alphaNew[iSolution], 3, 8) + " deg")

			iIteration += 1

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

def __TestOptimizeClimb__():
	lookupCl, lookupCd = TrackCommon.LoadAerodynamicData("./data/aerodynamicPerformance/Cl.csv",
														 "./data/aerodynamicPerformance/Cd.csv")
	OptimizeClimb(45000, 46500, 40000, 40, 0, 0, 0, 700*8.8, 35.0, -20, 80000, 0,
				  0.05, lookupCl, lookupCd)

__TestOptimizeClimb__()
