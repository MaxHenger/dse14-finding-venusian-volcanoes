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
import TrackBiasMap
import TrackStorage

def OptimizeDive(heightUpper, heightTarget, vHorInitial, vVerInitial,
                 longitude, latitude, W, S, vHorTarget, dt,
                 lookupCl, lookupCd, storeResults=True):
    atmosphere = Atmosphere.Atmosphere()

    # Retrieve ranges of angle of attack from the lookup tables
    numAlpha = 250
    alphaLimits = lookupCl.getPoints()
    alphaLimits = [alphaLimits[0][0], alphaLimits[0][-1]]

    # Settings for the optimization routine
    biasLimit = 0.1 # percent
    biasStep = 0.05 # percent
    biasWidth = 5000 # meters
    percentSpeedOfSound = 0.65

    alphaDotLimit = 0.5 # deg/s
    gammaDotLimit = 1.5 / 180.0 * np.pi # rad/s
    gammaLimit = np.pi / 2.0

    # Set initial values
    initialZonal = atmosphere.velocityZonal(heightUpper, latitude, longitude)[1]
    initialGamma = np.arctan2(-vVerInitial, initialZonal + vHorInitial)
    alpha = [0.0]
    vHor = [vHorInitial]
    vVer = [vVerInitial]
    height = [heightUpper]
    gamma = [initialGamma]
    time = [0.0]
    vLim = [0.0]

    # Set bias maps and associated variables
    biasBaseGamma = 0.975
    biasBaseVInf = 0.975
    biasBaseGammaDot = 0.975
    biasChooseGamma = 75.0 / 180.0 * np.pi # special variable: if vInf or gammaDot
        # exceed the allowed maxima but gamma is larger than this value, then
        # the gamma bias map will be adjusted, not the other maps

    biasHeightUpper = heightUpper + 0.2 * (heightUpper - heightTarget)
    biasHeightLower = heightTarget - 0.2 * (heightUpper - heightTarget)
    biasGamma = TrackBiasMap.BiasMap("gamma", biasHeightUpper, biasHeightLower, 1024, biasBaseGamma)
    biasVInf = TrackBiasMap.BiasMap("vInf", biasHeightUpper, biasHeightLower, 1024, biasBaseVInf)
    biasGammaDot = TrackBiasMap.BiasMap("gammaDot", biasHeightUpper, biasHeightLower, 1024, biasBaseGammaDot)

    # Start iterating
    solved = False

    while not solved:
        totalTime = 0.0

        alpha = [0.0]
        vHor = [vHorInitial]
        vVer = [vVerInitial]
        height = [heightUpper]
        gamma = [initialGamma]
        time = [0.0]
        vLim = [0.0]

        gammaOld = gamma[0]

        iIteration = 0
        solved = True

        while height[-1] > heightTarget:
            # Determine new valid range of angles of attack
            alphaNew = np.linspace(np.max([alphaLimits[0], alpha[-1] - dt * alphaDotLimit]),
                                   np.min([alphaLimits[1], alpha[-1] + dt * alphaDotLimit]), numAlpha)

            # Determine new flight variables
            vHorNew, vVerNew, gammaNew, hNew = TrackDive.Step(height[-1], alpha[-1],
                gamma[-1], vHor[-1], vVer[-1], longitude, latitude, W, S, alphaNew,
                dt, lookupCl, lookupCd, atmosphere, tol=1e-8, relax=0.8)

            totalTime = dt * (iIteration + 1)

            # Filter the valid solutions from the invalid ones. Keep track of
            # the reason why certain solutions fail in case all of them fail
            iValid = []
            vZonal = atmosphere.velocityZonal(hNew, latitude, longitude)[1]
            vInf = np.sqrt(np.power(vHorNew + vZonal, 2.0) + np.power(vVerNew, 2.0))
            vLimit = atmosphere.speedOfSound(hNew, latitude, longitude) * percentSpeedOfSound

            gammaOffenders = 0
            vInfOffenders = 0
            gammaDotOffenders = 0

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
                    # Sadly the initial gamma is quite oscillatory. Hence allow
                    # it to stabilize in the first few iterations
                    if abs((gammaNew[i] - gammaOld) / dt) > gammaDotLimit:
                        gammaDotOffenders += 1
                        isOffender = True

                if isOffender:
                    continue

                # This item is valid
                iValid.append(i)

            if len(iValid) == 0:
                # No valid solutions
                print("Did not find a solution at t =", totalTime)
                toShow = int(len(alphaNew) / 2)
                print(TrackCommon.StringPad("gamma  = ", gammaNew[toShow] * 180.0 / np.pi, 3, 10) + " deg")
                print(TrackCommon.StringPad("vInf   = ", vInf[toShow], 3, 10) + " m/s")
                print(TrackCommon.StringPad("vLimit = ", vLimit[toShow], 3, 10) + " m/s")
                print(TrackCommon.StringPad("vZonal = ", vZonal[toShow], 3, 10) + " m/s")
                print(TrackCommon.StringPad("vHor   = ", vHorNew[toShow], 3, 10) + " m/s")
                print(TrackCommon.StringPad("vVer   = ", vVerNew[toShow], 3, 10) + " m/s")
                print(TrackCommon.StringPad("gammaDot = ", abs((gammaNew[toShow] - gammaOld) * 180.0 / np.pi / dt), 5, 10) + " deg/s")
                print(TrackCommon.StringPad("gammaDot limit = ", gammaDotLimit * 180.0 / np.pi, 5, 10) + " deg/s")

                # Determine how to adjust the bias maps
                listOffenders = [gammaOffenders, vInfOffenders, gammaDotOffenders]
                iWorstOffender = np.argmax(listOffenders)

                if listOffenders[iWorstOffender] == 0:
                    print('No offenders, adjusting base biases:')
                    biasBaseGamma, biasBaseVInf, biasBaseGammaDot = \
                        TrackCommon.AdjustBiasMapCommonly([biasGamma, biasVInf, biasGammaDot],
                                                          biasStep, ['gamma', 'vInf', 'gammaDot'])
                else:
                    # For the reason behind the following indices, see the
                    # listOffenders variable declared above
                    print('average gamma:', np.average(abs(gammaNew)) * 180.0 / np.pi)
                    if iWorstOffender == 0 or np.average(abs(gammaNew)) * 180.0 / np.pi > biasChooseGamma:
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
                    else:
                        raise RuntimeError("Unrecognized offender index for bias map")

                if biasBaseGamma < biasLimit or biasBaseVInf < biasLimit or biasBaseGammaDot < biasLimit:
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
                # - base metric: maximizing vertical speed
                metric = ((vVerNew[i] + vLimit[i]) / vLimit[i])**2.0

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

                if metric < bestMetric:
                    bestMetric = metric
                    iSolution = i

            # Append the new values to the solution arrays
            height.append(hNew[iSolution])
            alpha.append(alphaNew[iSolution])
            gamma.append(gammaNew[iSolution])
            vHor.append(vHorNew[iSolution])
            vVer.append(vVerNew[iSolution])
            vLim.append(vLimit[iSolution])
            time.append(totalTime)

            gammaOld = gammaNew[iSolution]

            if iIteration % 30 == 0:
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

    axHor.plot(time, vHor, 'g')
    axHor.plot(timeFinal, vHorFinal, 'r')
    axHor.set_xlabel('time [s]')
    axHor.set_ylabel('horizontal speed [m/s]')
    axHor.grid(True)

    # Prepare Vinf
    vInf = np.zeros([len(vVer)])

    for i in range(0, len(vVer)):
        vZonal = atmosphere.velocityZonal(height[i], latitude, longitude)[1]
        vInf[i] = np.sqrt(np.power(vZonal + vHor[i], 2.0) + np.power(vVer[i], 2.0))

    VinfFinal = np.zeros([len(vVerFinal)])

    for i in range(0, len(vVerFinal)):
        vZonal = atmosphere.velocityZonal(heightFinal[i], latitude, longitude)[1]
        VinfFinal[i] = np.sqrt(np.power(vZonal + vHorFinal[i], 2.0) + np.power(vVerFinal[i], 2.0))

    axVinf.plot(time, vInf, 'g')
    axVinf.plot(timeFinal, VinfFinal, 'r')
    axVinf.plot(time, vLim, 'g--')
    axVinf.plot(timeFinal, vLimFinal, 'r--')
    axVinf.set_xlabel('time [s]')
    axVinf.set_ylabel('V_inf [m/s]')
    axVinf.grid(True)

    # Plot the bias maps
    fig = plt.figure()
    axBias = fig.add_subplot(111)
    lGamma, = axBias.plot(biasGamma.getAxis() / 1e3, biasGamma.getMap() * 1e2, 'r')
    lVInf, = axBias.plot(biasVInf.getAxis() / 1e3, biasVInf.getMap() * 1e2, 'g')
    lGammaDot, = axBias.plot(biasGammaDot.getAxis() / 1e3, biasVInf.getMap() * 1e2, 'b')
    axBias.legend([lGamma, lVInf, lGammaDot], ['gamma', 'vInf', 'gammaDot'])
    axBias.set_xlabel('h [km]')
    axBias.set_ylabel('bias [%]')
    axBias.grid(True)

    # Store results in a file
    if (storeResults):
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
        file.addVariable('biasGammaDot', biasGammaDot.getMap(), [biasGammaDot.getAxis()])
        file.addVariable('dt', dt)
        file.addVariable('vLim', vLim)
        file.save('dive' + str(heightUpper) + '-' + str(heightTarget) +
                  '.' + str(vHorInitial) + '-' + str(vHorTarget) + '.dat')

def PlotDive(filename):
    file = TrackStorage.DataStorage();
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
    dt = file.getVariable("dt").getValues()

    fig = plt.figure()
    axSpeed = fig.add_subplot(221)
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
    axHeight.plot(time, height, 'r')

    axHeight.set_xlabel(r'$t\;[s]$')
    axHeight.set_ylabel(r'$h\;[m]$')
    #axHeight.legend()
    axHeight.grid(True)

    # Plot the angles
    axAnglesLeft.plot(time, alpha, 'r', label=r'\alpha')
    for tick in axAnglesLeft.get_yticklabels():
        tick.set_color('r')

    axAnglesRight.plot(time, gamma * 180.0 / np.pi, 'g', label=r'\alpha')
    for tick in axAnglesRight.get_yticklabels():
        tick.set_color('g')

    axSpeed.set_xlabel(r'$t\;[s]$')
    axAnglesLeft.set_ylabel(r'$\alpha\;[\degree]$')
    axAnglesRight.set_ylabel(r'$\gamma\;[\degree]$')
    axAnglesLeft.grid(True)

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
    ax.set_xlabel(r'$h\;[km]$')
    ax.set_ylabel(r'$b\;[\%]$')
def __TestOptimizeDive__():
    lookupCl, lookupCd = TrackCommon.LoadAerodynamicData('./data/aerodynamicPerformance/Cl.csv',
                                                         './data/aerodynamicPerformance/Cd.csv')
    #OptimizeDive(55000, 35000, -20, -10, 0, 0, 700*8.8, 35.0, 0, 0.25, lookupCl, lookupCd)
    OptimizeDive(55000, 30000, -20, -10, 0, 0, 700*8.8, 35.0, 0, 0.25, lookupCl, lookupCd)
    #OptimizeDive(55000, 38000, -20, -10, 0, 0, 700*8.8, 35.0, 0, 0.25, lookupCl, lookupCd)
    #OptimizeDive(55000, 46000, -20, -10, 0, 0, 700*8.8, 35.0, 0, 0.25, lookupCl, lookupCd)
    #OptimizeDive(55000, 35000, -20, -10, 0, 0, 700*8.8, 35.0, 0, 0.25, lookupCl, lookupCd)
    #OptimizeDive(38000, 30000, 20, -10, 0, 0, 700*8.8, 35.0, 0, 0.25, lookupCl, lookupCd)
    #OptimizeDive(46000, 30000, 35, -10, 0, 0, 700*8.8, 35.0, 0, 0.25, lookupCl, lookupCd)

def __TestPlotDive__():
    PlotDive("dive55000-30000.-20-0.dat")

#__TestOptimizeDive__()
__TestPlotDive__()
