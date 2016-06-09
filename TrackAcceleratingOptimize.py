# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 23:03:22 2016

@author: MaxHenger
"""

import Atmosphere
import TrackAccelerating
import TrackCommon
import TrackBiasMap
import TrackAngleOfAttack

import numpy as np
import matplotlib.pyplot as plt

def OptimizeAccelerating(height, vHorInitial, PRequiredInitial, alphaInitial,
                         longitude, latitude, W, S, vHorFinal, inclination, dt,
                         PRequiredMin, PRequiredMax, lookupCl, lookupCd,
                         severity=0.0, plotResults=False, storeResults=True):
    # Construct lookup tables and interpolators
    atmosphere = Atmosphere.Atmosphere()
    lookupdCldAlpha = lookupCl.getDerivative()
    lookupdCddAlpha = lookupCd.getDerivative()
    g = 8.8

    # Set ranges in which the optimization routine can operate
    numPRequired = 125
    alphaLimits = lookupCl.getPoints()
    alphaLimits = [alphaLimits[0][0], alphaLimits[0][-1]]

    # Settings for the optimization routine
    biasLimit = 0.1
    biasStep = 0.15
    biasWidth = 5 #m/s
    percentSpeedOfSound = 0.75

    alphaDotLimit = 1.5 #degree/s
    PReqDotLimit = 250 #Watt/s
    updateCount = 35 # number of iterations to take before printing intermediate results

    steadySpeedValid = 50.0 # number of seconds to keep velocity steady to be valid
    steadySpeedRange = 1.0 # one-side of two-sided radius wherein speed should stay

    # Set initial values
    rho = TrackCommon.AdjustSeverity(atmosphere.density(height, latitude, longitude), severity)
    vZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(height, latitude, longitude), severity)
    initialVInf = vHorInitial + vZonal
    vLimit = atmosphere.speedOfSound(height, latitude, longitude) * percentSpeedOfSound

    # Retrieve expected final values
    alphaFinal, thrustFinal, valid = TrackAngleOfAttack.AngleOfAttackThrustSteady(
        W, S, 0.5 * rho * (vZonal + vHorFinal)**2.0, inclination, lookupCl,
        lookupdCldAlpha, lookupCd, lookupdCddAlpha)

    if valid == False:
        print(TrackCommon.StringPad("alpha final  = ", alphaFinal, 4, 8) + " deg")
        print(TrackCommon.StringPad("thrust final = ", thrustFinal, 4, 8) + " N")
        raise ValueError("Failed to obtain stable final thrust and angle of attack")

    PReqFinal = thrustFinal * (vZonal + vHorFinal)

    # Set up storage arrays
    alpha = [alphaInitial]
    vHor = [vHorInitial]
    power = [PRequiredInitial]
    time = [0.0]

    # Do some sanity checks
    if initialVInf > vLimit:
        raise ValueError("initial Vinf exceeds Vinf limit")

    if vHorFinal + vZonal > vLimit:
        raise ValueError("final vInf exceeds vInf limit")

    # Set bias maps and their associated values
    biasBaseDefaultAlpha = 0.975
    biasBaseDefaultAlphaDot = 0.975
    biasBaseDefaultVInf = 0.975

    biasBaseAlpha = biasBaseDefaultAlpha
    biasBaseAlphaDot = biasBaseDefaultAlphaDot
    biasBaseVInf = biasBaseDefaultVInf

    velocityBiasLower = vHorInitial - (vHorFinal - vHorInitial)
    velocityBiasUpper = vHorFinal + (vHorFinal - vHorInitial)

    biasAlpha = TrackBiasMap.BiasMap("alpha", velocityBiasLower, velocityBiasUpper, 1024, biasBaseDefaultAlpha)
    biasVInf = TrackBiasMap.BiasMap("vInf", velocityBiasLower, velocityBiasUpper, 1024, biasBaseDefaultVInf)
    biasAlphaDot = TrackBiasMap.BiasMap("alphaDot", velocityBiasLower, velocityBiasUpper, 1024, biasBaseDefaultAlphaDot)

    # Start iterating
    print(TrackCommon.StringHeader("Performing acceleration", 60))
    print(TrackCommon.StringPad(" > initial vHor: ", vHorInitial, 3, 10) + " m/s")
    print(TrackCommon.StringPad(" > initial vInf: ", initialVInf, 3, 10) + " m/s")
    print(TrackCommon.StringPad(" > final vHor:   ", vHorFinal, 3, 10) + " m/s")
    print(TrackCommon.StringPad(" > final vInf:   ", vHorFinal + vZonal, 3, 10) + " m/s")

    solved = False
    failed = False

    while not solved:
        totalTime = 0.0

        alpha = [alphaInitial]
        vHor = [vHorInitial]
        power = [PRequiredInitial]
        time = [0.0]

        alphaOld = alpha[0]

        iIteration = 0
        solved = True
        timeSteady = 0.0

        while timeSteady < steadySpeedValid:
            # Determine the new valid range of power values
            PReqNew = np.linspace(max(PRequiredMin, power[-1] - dt * PReqDotLimit),
                                  min(PRequiredMax, power[-1] + dt * PReqDotLimit), numPRequired)

            alphaNew, vHorNew, valid = TrackAccelerating.Step(vHor[-1], alpha[-1],
                longitude, latitude, power[-1], W, S, inclination, PReqNew, dt, lookupCl,
                lookupCd, lookupdCldAlpha, rho, g, vZonal)
            vInfNew = vHorNew + vZonal

            totalTime = dt * (iIteration + 1)
            alphaDot = (alphaNew - alphaOld) / dt

            # Filter the valid solutions from the invalid solutions. Keep track
            # of the offenders to possibly adjust the bias maps
            iValid = []

            alphaOffenders = 0
            vInfOffenders = 0
            alphaDotOffenders = 0

            for i in range(0, len(alphaNew)):
                # Determine the validity of this particular solution
                isOffender = False

                if alphaNew[i] < alphaLimits[0] or alphaNew[i] > alphaLimits[1]:
                    alphaOffenders += 1
                    isOffender = True

                if vInfNew[i] > vLimit:
                    vInfOffenders += 1
                    isOffender = True

                if abs(alphaDot[i]) > alphaDotLimit:
                    alphaDotOffenders += 1
                    isOffender = True

                if isOffender:
                    continue

                # This solution is a valid one, exciting!
                iValid.append(i)

            if len(iValid) == 0:
                # No valid solutions were found
                print('\n * Did not find a solution at:\n > ' +
                    TrackCommon.StringPad('t    = ', totalTime, 3, 10) + ' s\n')

                toShow = int(len(alphaNew) / 2)
                print(TrackCommon.StringPad(' > alpha           = ', alphaNew[toShow], 3, 10) + ' deg')
                print(TrackCommon.StringPad(' > alpha min limit = ', alphaLimits[0], 3, 10) + ' deg')
                print(TrackCommon.StringPad(' > alpha max limit = ', alphaLimits[1], 3, 10) + ' deg')
                print(TrackCommon.StringPad(' > alphaDot        = ', alphaDot[toShow], 3, 10) + ' deg/s')
                print(TrackCommon.StringPad(' > alphaDot limit  = ', alphaDotLimit, 3, 10) + ' deg/s')
                print(TrackCommon.StringPad(' > vHor            = ', vHorNew[toShow], 3, 10) + ' m/s')
                print(TrackCommon.StringPad(' > vInf            = ', vHorNew[toShow] + vZonal, 3, 10) + ' m/s')
                print(TrackCommon.StringPad(' > vInf limit      = ', vLimit, 3, 10) + ' m/s')

                # Determine how to adjust the bias maps
                listOffenders = [alphaOffenders, vInfOffenders, alphaDotOffenders]
                iWorstOffender = np.argmax(listOffenders)

                if listOffenders[iWorstOffender] == 0:
                    print('\n * No offenders, adjusting base biases')
                    biasBaseAlpha, biasBaseVInf, biasBaseAlphaDot = \
                        TrackCommon.AdjustBiasMapCommonly([biasAlpha, biasVInf,
                        biasAlphaDot], biasStep, ['alpha', 'vInf', 'alphaDot'])

                    print('')
                else:
                    # Adjust the bias maps locally
                    if iWorstOffender == 0:
                        # Adjust alpha bias map
                        biasBaseAlpha = TrackCommon.AdjustBiasMapIndividually(biasAlpha,
                            biasStep, vHor[-1], biasWidth, 'alpha')
                    elif iWorstOffender == 1:
                        # Adjust vInf bias map
                        biasBaseVInf = TrackCommon.AdjustBiasMapIndividually(biasVInf,
                            biasStep, vHor[-1], biasWidth, 'vInf')
                    elif iWorstOffender == 2:
                        # Adjust alphaDot bias map
                        biasBaseAlphaDot = TrackCommon.AdjustBiasMapIndividually(biasAlphaDot,
                            biasStep, vHor[-1], biasWidth, 'alphaDot')
                    else:
                        raise RuntimeError("Unrecognized offender index for bias map")

                    print('')

                if biasBaseAlpha < biasLimit or biasBaseVInf < biasLimit or \
                        biasBaseAlphaDot < biasLimit:
                    print("Failed to find a solution")
                    failed = True
                    break

                # Restart with a new bias
                solved = False
                break

            # If this code is reached then valid solutions exist. Filter the
            # best one using a least-squares method
            bestMetric = 1e19
            iSolution = 0

            for i in iValid:
                # Base metric: proximity to the final speed
                metric = ((vHorNew[i] - vHorFinal) / (vHorFinal - vHorInitial))**2.0
                metric += ((PReqNew[i] - PReqFinal) / (PRequiredMax))**2.0 * 0.00025

                # Influence of the bias maps
                # - proximity to the limiting angle of attack
                curBiasAlpha = biasAlpha(vHorNew[i])
                if alphaNew[i] < alphaLimits[0] * curBiasAlpha:
                    metric += ((alphaLimits[0] * curBiasAlpha - alphaNew[i]) /
                        (alphaLimits[0] * (1.0 - curBiasAlpha)))**2.0
                elif alphaNew[i] > alphaLimits[1] * curBiasAlpha:
                    metric += ((alphaNew[i] - curBiasAlpha * alphaLimits[1]) /
                        (alphaLimits[1] * (1.0 - curBiasAlpha)))**2.0

                # - proximity to the allowed freestream velocity
                curBiasVInf = biasVInf(vHorNew[i])
                if vInfNew[i] > vLimit * curBiasVInf:
                    metric += ((vInfNew[i] - vLimit * curBiasVInf) /
                        (vLimit * (1.0 - curBiasVInf)))**2.0

                # - proximity to the allowed change in angle of attack
                curBiasAlphaDot = biasAlphaDot(vHorNew[i])
                if abs(alphaDot[i]) > curBiasAlphaDot * alphaDotLimit:
                    metric += ((alphaDot[i] - curBiasAlphaDot * alphaDotLimit) /
                        (alphaDotLimit * (1.0 - curBiasAlphaDot)))**2.0

                if metric < bestMetric:
                    iSolution = i
                    bestMetric = metric

            # Append the new values to the solution arrays
            alpha.append(alphaNew[iSolution])
            vHor.append(vHorNew[iSolution])
            power.append(PReqNew[iSolution])
            time.append(totalTime)

            if abs(vHorNew[iSolution] - vHorFinal) < steadySpeedRange:
                timeSteady += dt
            else:
                timeSteady = 0.0

            alphaOld = alphaNew[iSolution]

            if iIteration % updateCount == 0:
                print(TrackCommon.StringPad("Solved at t = ", totalTime, 3, 8) +
                      TrackCommon.StringPad(" s, alpha = ", alphaNew[iSolution], 3, 5) +
                      TrackCommon.StringPad(" deg, vHor = ", vHorNew[iSolution], 2, 6) +
                      TrackCommon.StringPad(" m/s, PReq = ", PReqNew[iSolution] / 1e3, 3, 6) +
                      " kW")

            iIteration += 1

    # Plot results
    if plotResults == True:
        fig = plt.figure()
        axAlpha = fig.add_subplot(311)
        axVHor = fig.add_subplot(312)
        axPower = fig.add_subplot(313)

        axAlpha.plot(time, alpha, 'g')
        axAlpha.plot([time[0], time[-1]], [alphaLimits[0], alphaLimits[0]], 'k--')
        axAlpha.plot([time[0], time[-1]], [alphaLimits[-1], alphaLimits[-1]], 'k--')
        axAlpha.set_xlabel('time [s]')
        axAlpha.set_ylabel('alpha [deg]')
        axAlpha.grid(True)

        axVHor.plot(time, vHor, 'g')
        axVHor.plot([time[0], time[-1]], [vHorInitial, vHorInitial], 'k--')
        axVHor.plot([time[0], time[-1]], [vHorFinal, vHorFinal], 'k--')
        axVHor.set_xlabel('time [s]')
        axVHor.set_ylabel('vHor [m/s]')
        axVHor.grid(True)

        axPower.plot(time, np.asarray(power) / 1e3, 'g')
        axPower.plot([time[0], time[-1]], [PRequiredMin / 1e3, PRequiredMin / 1e3], 'k--')
        axPower.plot([time[0], time[-1]], [PRequiredMax / 1e3, PRequiredMax / 1e3], 'k--')
        axPower.set_xlabel('time [s]')
        axPower.set_ylabel('pReq [kW]')
        axPower.grid(True)

    return time, vHor, alpha, power

def TestAccelerating():
    lookupCl, lookupCd = TrackCommon.LoadAerodynamicData("./data/aerodynamicPerformance/Cl.csv",
                                                         "./data/aerodynamicPerformance/Cd.csv")
    OptimizeAccelerating(68e3, -30, 20e3, 9.0, 0, 0, 700*8.8, 35, 20, 0, 0.10,
        15e3, 32e3, lookupCl, lookupCd)

def TestDecelerating():
    lookupCl, lookupCd = TrackCommon.LoadAerodynamicData("./data/aerodynamicPerformance/Cl.csv",
                                                         "./data/aerodynamicPerformance/Cd.csv")
    OptimizeAccelerating(68e3, 20, 25e3, 2.3, 0, 0, 700*8.8, 35, -20, 0, 0.10,
        5e3, 32e3, lookupCl, lookupCd, plotResults=True)

def TestSpecialCase():
    lookupCl, lookupCd = TrackCommon.LoadAerodynamicData("./data/aerodynamicPerformance/Cl.csv",
                                                         "./data/aerodynamicPerformance/Cd.csv")
    OptimizeAccelerating(62e3, 6.976, 32e3, 0.847, 0, 0, 700*8.8, 35, 3.5, 0, 0.10,
        0e3, 32e3, lookupCl, lookupCd, 0, plotResults=True)

#TestAccelerating()
#TestDecelerating()
#TestSpecialCase()