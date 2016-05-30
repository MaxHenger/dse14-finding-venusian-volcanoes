# -*- coding: utf-8 -*-
"""
Created on Sun May 29 19:02:47 2016

@author: MaxHenger
"""

import matplotlib.pyplot as plt
import numpy as np

import Atmosphere
import TrackDive
import TrackCommon

def OptimizeDive(heightUpper, heightTarget, vHorInitial, vVerInitial,
                 longitude, latitude, W, S, vHorTarget, dt,
                 lookupCl, lookupCd):
    # Retrieve ranges of angle of attack from the lookup tables
    numAlpha = 351
    alphaNew = lookupCl.getPoints()
    alphaNew = np.linspace(alphaNew[0][0], alphaNew[0][-1], numAlpha)
    alphaNew = np.linspace(-2.5, 2.5, numAlpha)

    # Settings for the optimization routine
    biasLimit = 1.0
    biasStep = 0.15
    percentSpeedOfSound = 0.55
    weightVinf = 1.0
    weightGamma = 5.0
    alphaDotLimit = 10.0e5 * (alphaNew[1] - alphaNew[0]) / dt
    gammaDotLimit = 5.0e5
    gammaLimit = np.pi / 5.0

    # Set initial values
    alpha = [0.0]
    vHor = [vHorInitial]
    vVer = [vVerInitial]
    height = [heightUpper]
    gamma = [0.0]
    time = [0.0]
    vLim = [0.0]

    # Set global variables
    atmosphere = Atmosphere.Atmosphere()

    # Start iterating
    solved = False
    bias = 0.025

    while not solved:
        totalTime = 0.0

        alpha = [0.0]
        vHor = [vHorInitial]
        vVer = [vVerInitial]
        height = [heightUpper]
        gamma = [0.0]
        time = [0.0]
        vLim = [0.0]

        alphaOld = alpha[0]
        gammaOld = gamma[0]

        iIteration = 0
        solved = True

        while height[-1] > heightTarget:
            # Determine new flight variables
            vHorNew, vVerNew, gammaNew, hNew = TrackDive.Step(height[-1], alpha[-1],
                gamma[-1], vHor[-1], vVer[-1], longitude, latitude, W, S, alphaNew,
                dt, lookupCl, lookupCd, atmosphere, tol=1e-7, relax=0.4)

            totalTime = dt * (iIteration + 1)

            # Filter the valid solutions from the invalid ones
            iValid = []
            vZonal = atmosphere.velocityZonal(hNew, latitude, longitude)[1]
            vInf = np.sqrt(np.power(vHorNew + vZonal, 2.0) + np.power(vVerNew, 2.0))
            vLimit = atmosphere.speedOfSound(hNew, latitude, longitude) * \
                percentSpeedOfSound #* \
                #percentSpeedOfSound * (1 - bias)
            vVerMax = 0.0

            for i in range(0, len(alphaNew)):
                # Calculate local freestream velocity and the speed of sound
                #print('vVer =', vVerNew[i], ', gamma =', gammaNew[i] * 180.0 / np.pi,
                #      ', vInf =', np.sqrt(vInfSquared[i]), ', vLim =', vLimit[i])
                if gammaNew[i] > -gammaLimit and gammaNew[i] < gammaLimit and \
                        vInf[i] <= vLimit[i] and \
                        abs((alphaNew[i] - alphaOld) / dt) < alphaDotLimit and \
                        abs((gammaNew[i] - gammaOld) / dt) < gammaDotLimit:
                            
                    if vVerNew[i] > vVerMax:
                        vVerMax = vVerNew[i]

                    iValid.append(i)

            if len(iValid) == 0:
                # No valid solutions
                bias += biasStep

                print("Did not find a solution at t =", totalTime, ", bias = ", bias)
                print(TrackCommon.StringPad("gamma  = ", gammaNew[0] * 180.0 / np.pi, 8, 20))
                print(TrackCommon.StringPad("vInf   = ", vInf[0], 8, 20))
                print(TrackCommon.StringPad("vLimit = ", vLimit[0], 8, 20))
                print(TrackCommon.StringPad("vZonal = ", vZonal[0], 8, 20))
                print(TrackCommon.StringPad("vHor   = ", vHorNew[0], 8, 20))
                print(TrackCommon.StringPad("vVer   = ", vVerNew[0], 8, 20))

                if bias > biasLimit:
                    print("Failed to find a solution")
                    break

                # Restart with a new bias
                solved = False
                break

            # Take the valid solution where the vertical speed is maximized
            '''iSolution = iValid[0]
            vVerMax = vVerNew[iValid[0]]

            for i in range(1, len(iValid)):
                #print('valid vVer =', vVerNew[iValid[i]])
                if vVerNew[iValid[i]] < vVerMax:
                    vVerMax = vVerNew[iValid[i]]
                    iSolution = iValid[i]'''

            bestMetric = -1e9
            iSolution = 0

            for i in iValid:
                metric = (vVerNew[i] - vVerMax)**2.0

                metric /= (1 + (0.6*(alphaNew[i] - alphaOld)/dt)**2.0)

                #if alphaNew[i] * alphaOld < 0:
                #    metric /= 150.0e5

                metric /= (1 + (0.2 * ((vVerNew[i] - vVer[-1]) / vVerNew[-1]))**2.0)
                metric /= (1 + (0.2 * ((vHorNew[i] - vHor[-1]) / vHorNew[-1]))**2.0)

                if vInf[i] > vLimit[i] * (1.0 - bias):
                    metric *= ((vLimit[i] - vInf[i]) / (bias * vLimit[i]))**2.0
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

            print(TrackCommon.StringPad("Solved at t = ", totalTime, 3, 8) +
                  TrackCommon.StringPad(" s, h = ", hNew[iSolution], 0, 7) +
                  TrackCommon.StringPad(" m, Vver = ", vVerNew[iSolution], 2, 6) +
                  TrackCommon.StringPad(" m/s, gamma = ", gammaNew[iSolution] * 180.0 / np.pi, 3, 8) +
                  TrackCommon.StringPad(" deg, alpha = ", alphaNew[iSolution], 3, 8) + ' deg')

            iIteration += 1

        # If this position is reached then the algorithm iterated through all
        # height values

    alphaFinal = [0.0]
    vHorFinal = [vHorInitial]
    vVerFinal = [vVerInitial]
    heightFinal = [heightUpper]
    gammaFinal = [0.0]
    timeFinal = [0.0]
    vLimFinal = [0.0]

    if True:
        for iIteration in range(0, len(alpha)):
            alphaAverage = 0.0
            alphaMin = int(max(iIteration - 35, 0))
            alphaMax = int(min(iIteration + 35, len(alpha)))

            for i in range(alphaMin, alphaMax):
                alphaAverage += alpha[i]

            alphaAverage /= (alphaMax - alphaMin)

            vHorNew, vVerNew, gammaNew, hNew = TrackDive.Step(heightFinal[-1],
                alphaFinal[-1], gammaFinal[-1], vHorFinal[-1], vVerFinal[-1], longitude,
                latitude, W, S, np.asarray([alphaAverage]), dt, lookupCl, lookupCd, atmosphere)

            vLimFinal.append(atmosphere.speedOfSound(hNew[0], longitude, latitude) * percentSpeedOfSound)
            alphaFinal.append(alphaAverage)
            vHorFinal.append(vHorNew[0])
            vVerFinal.append(vVerNew[0])
            heightFinal.append(hNew[0])
            gammaFinal.append(gammaNew[0])
            timeFinal.append(time[iIteration])

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

    axGamma.plot(time, np.asarray(gamma) * 180.0 / np.pi, 'g')
    axGamma.plot(timeFinal, np.asarray(gammaFinal) * 180.0 / np.pi, 'r')
    axGamma.set_xlabel('time [s]')
    axGamma.set_ylabel('gamma [deg]')

    axHeight.plot(time, np.asarray(height) / 1e3, 'g')
    axHeight.plot(timeFinal, np.asarray(heightFinal) / 1e3, 'r')
    axHeight.set_xlabel('time [s]')
    axHeight.set_ylabel('height [km]')

    axVer.plot(time, vVer, 'g')
    axVer.plot(timeFinal, vVerFinal, 'r')
    axVer.set_xlabel('time [s]')
    axVer.set_ylabel('vertical speed [m/s]')

    axHor.plot(time, vHor, 'g')
    axHor.plot(timeFinal, vHorFinal, 'r')
    axHor.set_xlabel('time [s]')
    axHor.set_ylabel('horizontal speed [m/s]')

    # Prepare Vinf
    Vinf = np.zeros([len(vVer)])
    
    for i in range(0, len(vVer)):
        vZonal = atmosphere.velocityZonal(height[i], latitude, longitude)[1]
        Vinf[i] = np.sqrt(np.power(vZonal + vHor[i], 2.0) + np.power(vVer[i], 2.0))

    VinfFinal = np.zeros([len(vVerFinal)])
    
    for i in range(0, len(vVerFinal)):
        vZonal = atmosphere.velocityZonal(heightFinal[i], latitude, longitude)[1]
        VinfFinal[i] = np.sqrt(np.power(vZonal + vHorFinal[i], 2.0) + np.power(vVerFinal[i], 2.0))
        
    axVinf.plot(time, Vinf, 'g')
    axVinf.plot(timeFinal, VinfFinal, 'r')
    axVinf.plot(time, vLim, 'g--')
    axVinf.plot(timeFinal, vLimFinal, 'r--')
    axVinf.set_xlabel('time [s]')
    axVinf.set_ylabel('V_inf [m/s]')

def __TestOptimizeDive__():
    lookupCl, lookupCd = TrackCommon.LoadAerodynamicData('./data/aerodynamicPerformance/Cl.csv',
                                                         './data/aerodynamicPerformance/Cd.csv')
    OptimizeDive(55000, 35000, 0, -10, 0, 0, 700*8.8, 35.0, 0, 0.25, lookupCl, lookupCd)

__TestOptimizeDive__()
