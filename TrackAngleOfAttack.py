# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 18:49:23 2016

This file defines some common methods to determine the angle of attack
iteratively for several flight conditions. Every function has its own specific
assumptions to arrive at the generated values. Be sure to read them before using
them.

@author: MaxHenger
"""

import TrackCommon
import TrackIO

import numpy as np
import matplotlib.pyplot as plt # for quick debugging

# AngleOfAttackPowered will calculate the angle of attack in climbing
# accelerating flight with a given thrust such that the acceleration of the
# aircraft is pointed in the direction of the flight path angle.
def AngleOfAttackPowered(W, S, qInf, T, gamma, inclination, lookupCl, lookupdCldAlpha, Re, tol=1e-12):
	# Define some commonly used terms
	gammaSin = np.sin(gamma)
	gammaCos = np.cos(gamma)
	gammaSinCos = gammaSin * gammaCos

	toDeg = 180.0 / np.pi

	alphaMin = lookupCl.getPoints()[0][0] / toDeg
	alphaMax = lookupCl.getPoints()[0][-1] / toDeg

	# Define a lambda function for the normal root-finding angle of attack
	# function for readability of the code
	root = lambda W, S, qInf, T, alpha, gamma, inclination, lookupCl: \
		T * (np.cos(gamma + alpha + inclination) / gammaCos - \
			 np.sin(gamma + alpha + inclination) / gammaSin) - \
		qInf * S * lookupCl(alpha * toDeg, Re) / gammaSinCos + W / gammaCos

	# And a lambda function for the derivative of the root-finding function
	rootDeriv = lambda W, S, qInf, T, alpha, gamma, inclination, lookupdCldAlpha: \
		-T * (np.sin(gamma + alpha + inclination) / gammaCos +
			  np.cos(gamma + alpha + inclination) / gammaSin) - \
		qInf * S * lookupdCldAlpha(alpha * toDeg, Re) * toDeg / gammaSinCos

	# Perform iterations on an initial alpha = 0
	alphaOld = 0.0

	for i in range(0, 1000):
		alphaNew = alphaOld - root(W, S, qInf, T, alphaOld, gamma, inclination,
			lookupCl) / rootDeriv(W, S, qInf, T, alphaOld, gamma, inclination,
			lookupdCldAlpha)

		if alphaNew < alphaMin:
			return alphaMin * toDeg, False

		if alphaNew > alphaMax:
			return alphaMax * toDeg, False

		if abs(alphaNew - alphaOld) < tol:
			# Iterated towards stable solution
			return alphaNew * toDeg, True

		alphaOld = alphaNew

	# Failed to find a solution in the alotted iterations
	raise ValueError("Failed to iterate to a stable solution")
	#return 1337.0, False

# AngleOfAttackThrustClimbing calculates the angle of attack and the thrust
# necessary to facilitate climbing, nonaccelerating flight. As this algorithm
# can feature oscillations it will perform comparisons to a tolerance with
# a certain degree of allowed hysteresis
def AngleOfAttackThrustClimbing(W, S, qInf, gamma, inclination, lookupCl,
		lookupdCldAlpha, lookupCd, lookupdCddAlpha, Re, tol=1e-5):
	# Define the commonly used terms
	gammaSin = np.sin(gamma)
	gammaCos = np.cos(gamma)

	toDeg = 180.0 / np.pi

	alphaMin = lookupCl.getPoints()[0][0] / toDeg
	alphaMax = lookupCl.getPoints()[0][-1] / toDeg

	# Lambda functions (for code readability and reduction of redundancy)
	root = lambda W, S, qInf, alpha, gamma, inclination, lookupCl, lookupCd: \
		qInf * S * (
			lookupCl(alpha * toDeg, Re) * (gammaSin * np.tan(gamma + alpha + inclination) + gammaCos) +
			lookupCd(alpha * toDeg, Re) * (gammaCos * np.tan(gamma + alpha + inclination) - gammaSin)
		) - W

	rootDeriv = lambda W, S, qInf, alpha, gamma, inclination, lookupCl, \
			lookupdCldAlpha, lookupCd, lookupdCddAlpha: \
		qInf * S * (
			lookupdCldAlpha(alpha * toDeg, Re) * toDeg * (gammaSin * np.tan(gamma + alpha + inclination) + gammaCos) +
			lookupCl(alpha * toDeg, Re) * (gammaSin * 2.0 / (np.cos(2.0 * (gamma + alpha + inclination)) + 1.0)) +
			lookupdCddAlpha(alpha * toDeg, Re) * toDeg * (gammaCos * np.tan(gamma + alpha + inclination) + gammaSin) -
			lookupCd(alpha * toDeg, Re) * (gammaCos * 2.0 / (np.cos(2.0 * (gamma + alpha + inclination)) + 1.0))
		)

	thrust = lambda qInf, S, alpha, gamma, inclination, lookupCl, lookupCd: \
		qInf * S * (
			lookupCl(alpha * toDeg, Re) * gammaSin +
			lookupCd(alpha * toDeg, Re) * gammaCos
		) /  np.cos(gamma + alpha + inclination)

	# Perform iterations on an initial alpha = 0
	alphaOld = 0.0

	for i in range(0, 100):
		alphaNew = alphaOld - root(W, S, qInf, alphaOld, gamma, inclination,
			lookupCl, lookupCd) / rootDeriv(W, S, qInf, alphaOld, gamma,
			inclination, lookupCl, lookupdCldAlpha, lookupCd, lookupdCddAlpha)

		if alphaNew < 100 * alphaMin:
			return alphaMin * toDeg, 0, False

		if alphaNew > 100 * alphaMax:
			return alphaMax * toDeg, 0, False

		if abs(alphaNew - alphaOld) < tol:
			return alphaNew * toDeg, thrust(qInf, S, alphaNew, gamma,
				inclination, lookupCl, lookupCd), True

		#print('from', alphaOld, 'to', alphaNew)
		alphaOld = alphaNew

	# Failed to find a solution
	# plotAlpha = np.linspace(alphaMin, alphaMax, 1000)
	# plotRoot = np.zeros(plotAlpha.shape)
	#
	# for i in range(0, len(plotAlpha)):
	# 	plotRoot[i] = root(W, S, qInf, plotAlpha[i], gamma, inclination,
	# 		lookupCl, lookupCd)
	#
	# fig = plt.figure()
	# ax = fig.add_subplot(111)
	# ax.plot(plotAlpha * toDeg, plotRoot)
	# ax.set_xlabel('alpha [deg]')
	# ax.set_ylabel('root output [-]')
	# raise ValueError("Failed to iterate to a stable solution")
	return alphaMax * toDeg, 0, False

# AngleOfAttackThrustSteady performs the same calculations as
# AngleOfAttackThrustClimbing. This function, however, assumes that the flight
# path angle is 0. This simplifies the mathematics (and will cause a reduction
# in computation cost as a result).
def AngleOfAttackThrustSteady(W, S, qInf, inclination, lookupCl,
		lookupdCldAlpha, lookupCd, lookupdCddAlpha, Re, tol=1e-12):
	# Define the commonly used terms
	toDeg = 180.0 / np.pi

	alphaMin = lookupCl.getPoints()[0][0] / toDeg
	alphaMax = lookupCl.getPoints()[0][-1] / toDeg

	# Lambda functions
	root = lambda W, S, qInf, alpha, inclination, lookupCl, lookupCd: \
		qInf * S * (
			lookupCl(alpha * toDeg, Re) +
			lookupCd(alpha * toDeg, Re) * np.tan(alpha + inclination)
		) - W

	rootDeriv = lambda W, S, qInf, alpha, inclination, lookupdCldAlpha, lookupCd, lookupdCddAlpha: \
		qInf * S * (
			lookupdCldAlpha(alpha * toDeg, Re) * toDeg +
			lookupdCddAlpha(alpha * toDeg, Re) * toDeg * np.tan(alpha + inclination) +
			lookupCd(alpha * toDeg, Re) * 2.0 / (np.cos(2 * (alpha + inclination)) + 1)
		)

	thrust = lambda qInf, S, lookupCd, alpha, inclination: \
		qInf * S * lookupCd(alpha * toDeg, Re) / np.cos(alpha + inclination)

	# Perform iterations
	alphaOld = 0.0

	for i in range(0, 1000):
		alphaNew = alphaOld - root(W, S, qInf, alphaOld, inclination, lookupCl, lookupCd) / \
			rootDeriv(W, S, qInf, inclination, alphaOld, lookupdCldAlpha, lookupCd,
			lookupdCddAlpha)

		if alphaNew < alphaMin:
			return alphaMin * toDeg, 0, False

		if alphaNew > alphaMax:
			return alphaMax * toDeg, 0,  False

		if abs(alphaNew - alphaOld) < tol:
			return alphaNew * toDeg, thrust(qInf, S, lookupCd, alphaNew,
				inclination), True

		alphaOld = alphaNew

	# Failed to find a solution
	raise ValueError("Failed to iterate to a stable solution")
	#return 1337.0, 0, False

def AngleOfAttackSteady(W, S, qInf, lookupReverseCl, Re):
	res = lookupReverseCl(W / (qInf * S), Re)

	if len(res) != 1:
		return res[0], False

	return res[0], True

def __testAngleOfAttack__():
	# Do some simple testing
	lookupCl, lookupCd = TrackIO.LoadAerodynamicData("./data/aerodynamicPerformance/Cl.csv",
														 "./data/aerodynamicPerformance/Cd.csv")
	lookupdCldAlpha = lookupCl.getDerivative()
	lookupdCddAlpha = lookupCd.getDerivative()

	# Test the powered angle of attack
	print(' --- Testing angle of attack with specified thrust in climb')
	for qInf in [2, 3, 4]:
		aoa, found = AngleOfAttackPowered(1, 1, qInf, 1, 5 / 180.0 * np.pi, 0,
			lookupCl, lookupdCldAlpha)

		print(found, ', qInf =', qInf, ', aoa =', aoa)

		if found == False:
			raise ValueError("Failed to find plausible angle of attack:", aoa)

	# Do some testing on an impossible to attain angle of attack
	aoa, found = AngleOfAttackPowered(1, 1, 1e-10, 0, 5 / 180.0 * np.pi, 0,
		lookupCl, lookupdCldAlpha)

	if found == True:
		raise ValueError("Found implausible angle of attack:", aoa)

	# Test the angle of attack and thrust finding algorithm while climbing.
	print(' --- Testing angle of attack and thrust while climbing')
	for qInf in [2, 3, 4]:
		aoa, T, found = AngleOfAttackThrustClimbing(1, 1, qInf, 5 / 180.0 * np.pi, 0,
			lookupCl, lookupdCldAlpha, lookupCd, lookupdCddAlpha)

		print(found, ', qInf =', qInf, ', aoa =', aoa, ', T =', T)

		if found == False:
			raise ValueError("Failed to find plausible angle of attack:", aoa)

	aoa, T, found = AngleOfAttackThrustClimbing(1, 1, 1e-10, 5 / 180.0 * np.pi, 0,
		lookupCl, lookupdCldAlpha, lookupCd, lookupdCddAlpha)

	if found == True:
		raise ValueError("Found implausible angle of attack:", found)

	# Test the simplified steady angle of attack finding algorithm while climbing
	print(' --- Testing angle of attack and thrust while steady')
	for qInf in [2, 3, 4]:
		aoaRef, TRef, foundRef = AngleOfAttackThrustClimbing(1, 1, qInf, 0, 0,
			lookupCl, lookupdCldAlpha, lookupCd, lookupdCddAlpha)

		if foundRef == False:
			raise ValueError("Very weird, this used to work picoseconds ago")

		aoa, T, found = AngleOfAttackThrustSteady(1, 1, qInf, 0, lookupCl,
			lookupdCldAlpha, lookupCd, lookupdCddAlpha)

		print(found, ', qInf =', qInf, ', aoa =', aoa, ', T =', T)

		if found == False:
			raise ValueError("Failed to find plausible angle of attack:", aoa)

		if abs(TRef - T) > 1e-8:
			raise ValueError("Reference thrust", TRef, "differs too much from thrust", T)

		if abs(aoaRef - aoa) > 1e-8:
			raise ValueError("Reference aoa", aoaRef, "differs too much from aoa", aoa)

	aoa, T, found = AngleOfAttackThrustSteady(1, 1, 1e-10, 0, lookupCl,
		lookupdCldAlpha, lookupCd, lookupdCddAlpha)

	if found == True:
		raise ValueError("Found implausible angle of attack:", found)

#__testAngleOfAttack__()
