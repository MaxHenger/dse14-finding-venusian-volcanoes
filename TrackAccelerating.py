# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 16:32:43 2016

@author: MaxHenger
"""

import Atmosphere
import TrackLookup
import TrackCommon
import TrackIO

import matplotlib.pyplot as plt
import numpy as np

def AngleOfAttack(W, S, vInf, qInf, inclination, PNew, alphaMin, alphaMax,
		lookupCl, lookupdCldAlpha, tol=1e-8):
	# Define lambda functions
	toDeg = 180.0 / np.pi

	root = lambda W, S, vInf, qInf, alpha, inclination, P, lookupCl: \
		qInf * S * lookupCl(alpha * toDeg) +  P / vInf * np.sin(alpha + inclination) - W

	deriv = lambda S, vInf, qInf, alpha, inclination, P, lookupdCldAlpha: \
		qInf * S * lookupdCldAlpha(alpha * toDeg) * toDeg + P / vInf * np.cos(alpha + inclination)

	# Keep track of the non-converged angles of attack
	alphaOld = 0.0

	# Iterate until converged
	for i in range(0, 1000):
		alphaNew = alphaOld - root(W, S, vInf, qInf, alphaOld, inclination, PNew,
			lookupCl) / deriv(S, vInf, qInf, alphaOld, inclination, PNew,
			lookupdCldAlpha)

		if alphaNew < alphaMin:
			return alphaMin, False

		if alphaNew > alphaMax:
			return alphaMax, False

		if abs(alphaNew - alphaOld) < tol:
			return alphaNew, True

		#print('from', round(alphaOld * toDeg, 9), 'to', round(alphaNew * toDeg, 9))

		alphaOld = alphaNew

	# All points did not converge
	raise ValueError("TrackAccelerating.AngleOfAttack did not converge")

def Step(vHorCur, alphaCur, longitudeCur, latitudeCur, PCur,
		 W, S, inclination, PNew, dt, lookupCl, lookupCd, lookupdCldAlpha, rho,
		 g, vZonal, tol=1e-8, relax=0.8):
	# Some often used-variables
	toDeg = 180.0 / np.pi
	vInf1 = vZonal + vHorCur
	qInf1 = 0.5 * rho * vInf1**2.0
	term1 = PCur / vInf1 * np.cos(alphaCur / toDeg + inclination) - \
		qInf1 * S * lookupCd(alphaCur)

	clPoints = lookupCl.getPoints()
	alphaLimits = [clPoints[0][0] / toDeg, clPoints[0][-1] / toDeg]

	# Get the initial angles of attack and the initial new velocity
	alphaNew = []
	vHorNew = []
	validNew = []

	for i in range(0, len(PNew)):
		# Calculate initial values
		alpha, v = AngleOfAttack(W, S, vInf1, qInf1, inclination,
			PNew[i], alphaLimits[0], alphaLimits[1], lookupCl, lookupdCldAlpha)

		vHorOld = vHorCur + 0.5 * g / W * dt * ( term1
			+ PNew[i] / vInf1 * np.cos(alpha + inclination)
			- qInf1 * S * lookupCd(alpha * toDeg)
		)

		# Start iterating
		converged = False

		for j in range(0, 1000):
			# Calculate the new atmospheric properties and angle of attack
			vInf2 = vZonal + vHorCur
			qInf2 = 0.5 * rho * vInf2**2.0

			alpha, v = AngleOfAttack(W, S, vInf2, qInf2, inclination, PNew[i],
				alphaLimits[0], alphaLimits[1], lookupCl, lookupdCldAlpha)

			if not v:
				# Angle of attack solution was invalid
				break

			# Calculate the new horizontal velocity
			vHor = vHorCur + 0.5 * g / W * dt * ( term1
				+ PNew[i] / vInf2 * np.cos(alpha + inclination)
				- qInf2 * S * lookupCd(alpha * toDeg)
			)

			if abs((vHor - vHorOld) / vHorOld) < tol:
				alphaNew.append(alpha * toDeg)
				vHorNew.append(vHor)
				validNew.append(True)
				converged = True
				break

			vHorOld = vHor

		if not converged:
			alphaNew.append(alphaLimits[1] * toDeg)
			vHorNew.append(0)
			validNew.append(False)

	return np.asarray(alphaNew), np.asarray(vHorNew), np.asarray(validNew)

def TestStep(vHor, alpha, longitude, latitude, PCur, W, S, inclination, dt, h):
	PNew = np.linspace(0.0, 2.0 * PCur, 10)

	atm = Atmosphere.Atmosphere()
	rho = atm.density(h, latitude, longitude)[1]
	zonal = atm.velocityZonal(h, latitude, longitude)[1]

	lookupCl, lookupCd = TrackIO.LoadAerodynamicData("./data/aerodynamicPerformance/Cl.csv",
														 "./data/aerodynamicPerformance/Cd.csv")

	newAlpha, newVHor, valid = Step(vHor, alpha, longitude, latitude, PCur, W, S,
		inclination, PNew, dt, lookupCl, lookupCd, lookupCl.getDerivative(),
		rho, 8.8, zonal)

	print(TrackCommon.StringHeader('Lengths', 40))
	print('length PNew:    ', len(PNew))
	print('length newAlpha:', len(newAlpha))
	print('length newVHor: ', len(newVHor))
	print('length valid:   ', len(valid))

	print(TrackCommon.StringHeader('Output', 40))
	for i in range(0, len(PNew)):
		print(TrackCommon.StringPad("Power: ", PNew[i] / 1e3, 3, 8) + " kW, " +
			  TrackCommon.StringPad("alpha: ", newAlpha[i], 3, 8) + " deg, " +
			  TrackCommon.StringPad("vHor: ", newVHor[i], 3, 8) + " m/s, " +
			  "valid:", valid[i])

	fig = plt.figure()
	ax = fig.add_subplot(121)
	ax.plot(PNew, newAlpha, 'g')
	ax.plot([PNew[0], PNew[-1]], [alpha, alpha], 'k--')
	ax.set_xlabel('PNew [W]')
	ax.set_ylabel('alpha [deg]')

	ax = fig.add_subplot(122)
	ax.plot(PNew, newVHor, 'g')
	ax.plot([PNew[0], PNew[-1]], [vHor, vHor], 'k--')
	ax.set_xlabel("PNew [W]")
	ax.set_ylabel("vHor [m/s]")


#TestStep(20, 3, 0, 0, 32e3, 700*8.8, 35, 0, 1.25, 68e3)
