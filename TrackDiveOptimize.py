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
                 longitude, latitude, W, S, vHorTarget, vVerTarget, dt,
                 lookupCl, lookupCd, storeResults=True):
    # Variables ONLY used for debugging. All pieces of code referencing them
    # are prefixed with the 'DEBUG' term
    plotAndQuit = False
    numPlotted = 0
    maxPlotted = 15
    plotAndQuitAxes = []
    plotAndQuitCmap = None

    # Create an atmosphere object
    atmosphere = Atmosphere.Atmosphere()

    # Retrieve ranges of angle of attack from the lookup tables
    numAlpha = 300
    alphaLimits = lookupCl.getPoints()
    alphaLimits = [alphaLimits[0][0], alphaLimits[0][-1]]

    # Settings for the optimization routine
    biasLimit = 0.1 # percent
    biasStep = 0.15 # percent
    biasWidth = 5000 # meters
    percentSpeedOfSound = 0.65

    alphaDotLimit = 0.5 # deg/s
    gammaDotLimit = 1.5 / 180.0 * np.pi # rad/s
    gammaLimit = np.pi / 2.0 # maximum negative and positive diving angle

    flareFailHeight = heightUpper - (heightUpper - heightTarget) * 0.1
    flareGammaValid = 2.0 / 180.0 * np.pi # one side of a two-sided range in which the flare angle is acceptable
    flareHeightValid = 250 # one side of a two-sided range in which the final height is acceptable
    updateCount = 35 # number of iterations before printing an update statement

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
    defaultBiasBaseGamma = 0.975
    defaultBiasBaseVInf = 0.975
    defaultBiasBaseGammaDot = 0.75
    biasBaseGamma = defaultBiasBaseGamma
    biasBaseVInf = defaultBiasBaseVInf
    biasBaseGammaDot = defaultBiasBaseGammaDot
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

    print(TrackCommon.StringHeader("Optimizing diving", 60))

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
            gammaDot = (gammaNew - gammaOld) / dt

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
                    if abs(gammaDot[i]) > gammaDotLimit:
                        gammaDotOffenders += 1
                        isOffender = True

                if isOffender:
                    continue

                # This item is valid
                iValid.append(i)

            if len(iValid) == 0:
                # No valid solutions
                print('\n * Did not find a solution at:\n > ' +
                    TrackCommon.StringPad("t = ", totalTime, 3, 10) + ' s\n > ' +
                    TrackCommon.StringPad("h = ", hNew[-1], 3, 10) + ' m\n')

                #toShow = int(len(alphaNew) / 2)
                #print(TrackCommon.StringPad("gamma  = ", gammaNew[toShow] * 180.0 / np.pi, 3, 10) + " deg")
                #print(TrackCommon.StringPad("vInf   = ", vInf[toShow], 3, 10) + " m/s")
                #print(TrackCommon.StringPad("vLimit = ", vLimit[toShow], 3, 10) + " m/s")
                #print(TrackCommon.StringPad("vZonal = ", vZonal[toShow], 3, 10) + " m/s")
                #print(TrackCommon.StringPad("vHor   = ", vHorNew[toShow], 3, 10) + " m/s")
                #print(TrackCommon.StringPad("vVer   = ", vVerNew[toShow], 3, 10) + " m/s")
                #print(TrackCommon.StringPad("gammaDot = ", abs(gammaDot[toShow]) * 180.0 / np.pi), 5, 10) + " deg/s")
                #print(TrackCommon.StringPad("gammaDot limit = ", gammaDotLimit * 180.0 / np.pi, 5, 10) + " deg/s")

                # DEBUG: If 'plotAndQuit' is set to true, plot the first couple
                # of failing solutions and stop when 'numPlotted' equals
                # 'maxPlotted'
                if plotAndQuit == True:
                    if len(plotAndQuitAxes) == 0:
                        fig = plt.figure()
                        plotAndQuitAxes.append(fig.add_subplot(131))
                        plotAndQuitAxes.append(fig.add_subplot(132))
                        plotAndQuitAxes.append(fig.add_subplot(133))
                        yLabel = ['gamma [deg]', 'abs gammaDot [deg/s]', 'v [m/s]']

                        for iPlot in range(0, len(plotAndQuitAxes)):
                            plotAndQuitAxes[iPlot].set_xlabel('alpha [deg]')
                            plotAndQuitAxes[iPlot].set_ylabel(yLabel[iPlot])
                            plotAndQuitAxes[iPlot].grid(True)

                        plotAndQuitCmap = plt.get_cmap('jet')

                    lineColor = plotAndQuitCmap((numPlotted + 1) / maxPlotted)
                    lineLabel = 'attempt ' + str(numPlotted + 1) + " it " + str(iIteration)
                    plotAndQuitAxes[0].plot(alphaNew, gammaNew * 180.0 / np.pi, color=lineColor, label=lineLabel)
                    plotAndQuitAxes[1].plot(alphaNew, abs((gammaNew - gammaOld) * 180.0 / np.pi / dt), color=lineColor, label=lineLabel)
                    plotAndQuitAxes[2].plot(alphaNew, vInf, color=lineColor, label=lineLabel)

                    numPlotted += 1

                    if numPlotted == maxPlotted:
                        plotAndQuitAxes[0].legend()
                        plotAndQuitAxes[1].legend()
                        plotAndQuitAxes[2].legend()
                        plotAndQuitAxes[1].plot([alphaNew[0], alphaNew[-1]], [gammaDotLimit, gammaDotLimit], 'k--')
                        plotAndQuitAxes[2].plot([alphaNew[0], alphaNew[-1]], [vLimit, vLimit], 'k--')

                        # Stop running
                        return

                # Determine how to adjust the bias maps
                listOffenders = [gammaOffenders, vInfOffenders, gammaDotOffenders]
                iWorstOffender = np.argmax(listOffenders)

                if listOffenders[iWorstOffender] == 0:
                    print('\n * No offenders, adjusting base biases:')
                    biasBaseGamma, biasBaseVInf, biasBaseGammaDot = \
                        TrackCommon.AdjustBiasMapCommonly([biasGamma, biasVInf, biasGammaDot],
                                                          biasStep, ['gamma', 'vInf', 'gammaDot'])
                    print('')
                else:
                    # For the reason behind the following indices, see the
                    # listOffenders variable declared above
                    print('\n * Adjusting bias, average gamma:', np.average(abs(gammaNew)) * 180.0 / np.pi)

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
                    print('')

                if biasBaseGamma < biasLimit or biasBaseVInf < biasLimit or biasBaseGammaDot < biasLimit:
                    print("Failed to find a solution")
                    return
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
                curBiasGammaDot = biasGammaDot(hNew[i])
                if gammaDot[i] > curBiasGammaDot * gammaDotLimit:
                    metric += ((gammaDot[i] - gammaDotLimit * curBiasGammaDot) /
                        (gammaDotLimit * (1.0 - curBiasGammaDot)))**2.0

                # - influence of freestream velocity's proximity to the limit
                curBiasVInf = biasVInf(hNew[i])
                if vInf[i] > vLimit[i] * curBiasVInf:
                    metric += ((vInf[i] - vLimit[i] * curBiasVInf) /
                        (vLimit[i] * (1.0 - curBiasVInf)))**2.0

                # - influence of flight path angle's proximity to the limit
                curBiasGamma = biasGamma(hNew[i])
                if abs(gammaNew[i]) > curBiasGamma * gammaLimit:
                    metric += ((abs(gammaNew[i]) - curBiasGamma * gammaLimit) /
                        (gammaLimit * (1.0 - curBiasGamma)))**2.0

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

            if iIteration % updateCount == 0:
                print(TrackCommon.StringPad("Solved at t = ", totalTime, 3, 8) +
                      TrackCommon.StringPad(" s, h = ", hNew[iSolution], 0, 7) +
                      TrackCommon.StringPad(" m, Vver = ", vVerNew[iSolution], 2, 6) +
                      TrackCommon.StringPad(" m/s, gamma = ", gammaNew[iSolution] * 180.0 / np.pi, 3, 8) +
                      TrackCommon.StringPad(" deg, alpha = ", alphaNew[iSolution], 3, 8) + ' deg')

            iIteration += 1

        # If this position is reached then the algorithm iterated through all
        # height values

    # Start the phase where flaring is performed. Start flaring halfway in the
    # dive. If somehow the starting height for flaring comes within a certain
    # percentage of the upper height then the dive is considered impossible
    print(TrackCommon.StringHeader("Optimizing Flaring", 60))

    # ascertain final freestream velocity and flight path angle
    vZonalFinal = atmosphere.velocityZonal(heightTarget, latitude, longitude)[1]
    vInfFinal = np.sqrt(np.power(vHorTarget + vZonalFinal, 2.0) + np.power(vVerTarget, 2.0))
    gammaFinal = np.arctan2(-vVerTarget, vZonalFinal + vHorTarget)

    # setup heights for flaring iterations
    iterativeUpperHeight = heightUpper
    iterativeLowerHeight = heightTarget
    iterativeCenterHeight = (heightUpper + heightTarget) / 2.0

    # setup bias map base values for flaring
    biasBaseFlareGamma = defaultBiasBaseGamma
    biasBaseFlareVInf = defaultBiasBaseVInf
    biasBaseFlareGammaDot = defaultBiasBaseGammaDot

    biasFlareHeightUpper = iterativeCenterHeight + 0.2 * (iterativeUpperHeight - iterativeLowerHeight)
    biasFlareHeightLower = iterativeCenterHeight - 0.2 * (iterativeUpperHeight - iterativeLowerHeight)

    biasFlareGamma = TrackBiasMap.BiasMap("gamma", biasFlareHeightUpper, biasFlareHeightLower, 1024, biasBaseFlareGamma)
    biasFlareVInf = TrackBiasMap.BiasMap("vInf", biasFlareHeightUpper, biasFlareHeightLower, 1024, biasBaseFlareVInf)
    biasFlareGammaDot = TrackBiasMap.BiasMap("gammaDot", biasFlareHeightUpper, biasFlareHeightLower, 1024, biasBaseFlareGammaDot)

    print(TrackCommon.StringPad(" * Final height = ", heightTarget, 1, 10) + ' km')
    print(TrackCommon.StringPad(" * Final gamma  = ", gammaFinal * 180.0 / np.pi, 3, 10) + ' deg')
    print(TrackCommon.StringPad(" * Target vInf  = ", vInfFinal, 2, 10) + " m/s")

    while iterativeCenterHeight < flareFailHeight:
        # Find to with iteration step this height corresponds
        iStartIteration = 0

        for i in range(0, len(height)):
            if height[i] < iterativeCenterHeight:
                if i == 0:
                    raise ValueError("Initial height is smaller than center height. " +
                        "This should be impossible!")

                iStartIteration = i - 1
                break

        print(TrackCommon.StringPad("\n * Begin flare =  ", height[iStartIteration], 1, 10) + " m\n")

        # Setup initial values for flaring
        alphaFlare = [alpha[iStartIteration]]
        vHorFlare = [vHor[iStartIteration]]
        vVerFlare = [vVer[iStartIteration]]
        heightFlare = [height[iStartIteration]]
        gammaFlare = [gamma[iStartIteration]]
        timeFlare = [time[iStartIteration]]
        vLimFlare = [vLim[iStartIteration]]

        gammaOld = gammaFlare[0]

        iIteration = 0
        solved = True

        while heightFlare[-1] > heightTarget - flareHeightValid and \
                abs(gammaFlare[-1] - gammaFinal) > flareGammaValid:
            # Determine new valid range of angles of attack
            alphaNew = np.linspace(np.max([alphaLimits[0], alphaFlare[-1] - dt * alphaDotLimit]),
                                   np.min([alphaLimits[1], alphaFlare[-1] + dt * alphaDotLimit]), numAlpha)

            # Determine new flight variables
            vHorNew, vVerNew, gammaNew, hNew = TrackDive.Step(heightFlare[-1],
                alphaFlare[-1], gammaFlare[-1], vHorFlare[-1], vVerFlare[-1],
                longitude, latitude, W, S, alphaNew, dt, lookupCl, lookupCd,
                atmosphere, tol=1e-8, relax=0.8)

            totalTime = timeFlare[0] + dt * (iIteration + 1)

            # Filter the valid solutions from the invalid ones. Keep track of
            # why certain solutions fail
            iValid = []
            vZonal = atmosphere.velocityZonal(hNew, latitude, longitude)[1]
            vInf = np.sqrt(np.power(vHorNew + vZonal, 2.0) + np.power(vVerNew, 2.0))
            vLimit = atmosphere.speedOfSound(hNew, latitude, longitude) * percentSpeedOfSound
            gammaDot = (gammaNew - gammaOld) / dt

            gammaOffenders = 0
            vInfOffenders = 0
            gammaDotOffenders = 0

            for i in range(0, len(alphaNew)):
                isOffender = False

                if gammaNew[i] < -gammaLimit or gammaNew[i] > gammaLimit:
                    gammaOffenders += 1
                    isOffender = True

                if vInf[i] > vLimit[i]:
                    vInfOffenders += 1
                    isOffender = True

                if abs(gammaDot[i]) > gammaDotLimit:
                    gammaDotOffenders += 1
                    isOffender = True

                if isOffender:
                    continue

                # This item is valid
                iValid.append(i)

            if len(iValid) == 0:
                # No valid solutions
                print('\n * Did not find a solution at:\n > ' +
                    TrackCommon.StringPad("t = ", totalTime, 3, 10) + ' s\n > ' +
                    TrackCommon.StringPad("h = ", hNew[-1], 3, 10) + ' m\n')

                #toShow = int(len(alphaNew) / 2)
                #print(TrackCommon.StringPad("gamma  = ", gammaNew[toShow] * 180.0 / np.pi, 3, 10) + " deg")
                #print(TrackCommon.StringPad("vInf   = ", vInf[toShow], 3, 10) + " m/s")
                #print(TrackCommon.StringPad("vLimit = ", vLimit[toShow], 3, 10) + " m/s")
                #print(TrackCommon.StringPad("vZonal = ", vZonal[toShow], 3, 10) + " m/s")
                #print(TrackCommon.StringPad("vHor   = ", vHorNew[toShow], 3, 10) + " m/s")
                #print(TrackCommon.StringPad("vVer   = ", vVerNew[toShow], 3, 10) + " m/s")
                #print(TrackCommon.StringPad("gammaDot = ", abs(gammaDot[toShow]) * 180.0 / np.pi, 5, 10) + " deg/s")
                #print(TrackCommon.StringPad("gammaDot limit = ", gammaDotLimit * 180.0 / np.pi, 5, 10) + " deg/s")

                # Determine how to adjust the bias maps
                listOffenders = [gammaOffenders, vInfOffenders, gammaDotOffenders]
                iWorstOffender = np.argmax(listOffenders)

                if listOffenders[iWorstOffender] == 0:
                    print('\n * No offenders, adjusting base biases:')
                    biasBaseFlareGamma, biasBaseFlareVInf, biasBaseFlareGammaDot = \
                        TrackCommon.AdjustBiasMapCommonly([biasBaseFlareGamma,
                        biasBaseFlareVInf, biasBaseFlareGammaDot], biasStep,
                        ['gamma', 'vInf', 'gammaDot'])
                    print('')
                else:
                    # Figure out how to adjust the bias map
                    print(TrackCommon.StringPad('\n * Adjusting bias, average gamma: ',
                        np.average(abs(gammaNew)) * 180.0 / np.pi, 2, 10))

                    if iWorstOffender == 0:
                        biasBaseFlareGamma = TrackCommon.AdjustBiasMapIndividually(
                            biasFlareGamma, biasStep, heightFlare[-1], biasWidth, 'gamma')
                    elif iWorstOffender == 1:
                        biasBaseFlareVInf = TrackCommon.AdjustBiasMapIndividually(
                            biasFlareVInf,  biasStep, heightFlare[-1], biasWidth, 'vInf')
                    elif iWorstOffender == 2:
                        biasBaseFlareGammaDot = TrackCommon.AdjustBiasMapIndividually(
                            biasFlareGammaDot, biasStep, heightFlare[-1], biasWidth, 'gammaDot')
                    else:
                        raise RuntimeError("Unrecognized flare offender index for bias map")

                    print('')

                if biasBaseFlareGamma < biasLimit or biasBaseFlareVInf < biasLimit or \
                        biasBaseFlareGammaDot < biasLimit:
                    print("Failed to find a flaring solution")
                    return

                # Restart with a new bias
                solved = False
                break

            # Within the valid solutions seek the solution that maximizes a
            # metric defined to get the aircraft to the intended height at the
            # intended velocity and flight path angle
            bestMetric = 1e19
            iSolution = 0

            for i in iValid:
                # Base metric contributions
                # - closing in on the desired speed
                # NOTE: Experimental addition is the linear scaling based on the
                # distance to the final intended altitude
                metric = ((vInf[i] - vInfFinal) / vLimit[i])**2.0
                metric *= 2.5 * (1 - (hNew[i] - heightTarget) / (heightFlare[0] - heightTarget))

                # - closing in on the desired flight path angle
                metric += ((gammaNew[i] - gammaFinal) / gammaLimit)**2.0

                # Modifying contributions to steer away from limiting regions
                # - influence of dgamma/dt
                curBiasGammaDot = biasFlareGammaDot(hNew[i])
                if gammaDot[i] > curBiasGammaDot * gammaDotLimit:
                    metric += ((gammaDot[i] - gammaDotLimit * curBiasGammaDot) /
                        (gammaDotLimit * (1.0 - curBiasGammaDot)))**2.0

                # influence of the freestream velocity
                curBiasVInf = biasFlareVInf(hNew[i])
                if vInf[i] > vLimit[i] * curBiasVInf:
                    metric += ((vInf[i] - vLimit[i] * curBiasVInf) /
                        (vLimit[i] * (1.0 - curBiasVInf)))**2.0

                # influence of the flight path angle
                curBiasGamma = biasFlareGamma(hNew[i])
                if abs(gammaNew[i]) > curBiasGamma * gammaLimit:
                    metric += ((abs(gammaNew[i]) - curBiasGamma * gammaLimit) /
                        (gammaLimit * (1.0 - curBiasGamma)))**2.0

                if metric < bestMetric:
                    bestMetric = metric
                    iSolution = i

            # Append the new values to the solution arrays
            heightFlare.append(hNew[iSolution])
            alphaFlare.append(alphaNew[iSolution])
            gammaFlare.append(gammaNew[iSolution])
            vHorFlare.append(vHorNew[iSolution])
            vVerFlare.append(vVerNew[iSolution])
            vLimFlare.append(vLimit[iSolution])
            timeFlare.append(totalTime)

            # Check if the current solution adheres to both requirements
            if (abs(heightFlare[-1] - heightTarget) < flareHeightValid and
                    abs(gammaFlare[-1] - gammaFinal) < flareGammaValid):
                # Found a solution!
                print("\n * Solution found at:\n > " +
                    TrackCommon.StringPad("t     = ", totalTime, 3, 10) + ' s\n > ' +
                    TrackCommon.StringPad("h     = ", heightFlare[-1], 1, 10) + ' m\n > ' +
                    TrackCommon.StringPad("gamma = ", gammaFlare[-1] * 180.0 / np.pi, 2, 10) + " deg\n")
                break

            gammaOld = gammaNew[iSolution]

            if iIteration % updateCount == 0:
                print(TrackCommon.StringPad("Solved at t = ", totalTime, 3, 8) +
                      TrackCommon.StringPad(" s, h = ", hNew[iSolution], 0, 7) +
                      TrackCommon.StringPad(" m, Vver = ", vVerNew[iSolution], 2, 6) +
                      TrackCommon.StringPad(" m/s, gamma = ", gammaNew[iSolution] * 180.0 / np.pi, 3, 8) +
                      TrackCommon.StringPad(" deg, alpha = ", alphaNew[iSolution], 3, 8) + ' deg')

            iIteration += 1

        if solved == True:
            if abs(gammaFlare[-1] - gammaFinal) > flareGammaValid:
                # Flare angle differs by too much, but we've reached the intended
                # height.
                if gammaFlare[-1] < gammaFinal:
                    # Flight path angle is less steep then the final flight path
                    # angle, I suspect this will never happen
                    raise RuntimeError("Flight path angle was unexpectedly less "
                        + "steep when the target height was reached")
                else:
                    # Flight path angle is too steep, initiate the flare earlier
                    iterativeLowerHeight = iterativeCenterHeight
                    iterativeCenterHeight = (iterativeLowerHeight + iterativeUpperHeight) / 2.0
                    print('\n * Flight path angle too steep, initiating flare at',
                        round(iterativeCenterHeight / 1e3, 3), 'km\n')
            elif abs(heightFlare[-1] - heightTarget) > flareHeightValid:
                # Intended flight path angle is reached but the intended height
                # is not yet reached, the flaring can be initiated later
                iterativeUpperHeight = iterativeCenterHeight
                iterativeCenterHeight = (iterativeLowerHeight + iterativeUpperHeight) / 2.0
                print('\n * Height at which flare end is too high, initiating flare at',
                    round(iterativeCenterHeight / 1e3, 3), 'km\n')
            else:
                # Strip the old solutions of their dive that is not turned into
                # a flaring movement
                alpha = alpha[:iStartIteration]
                vHor = vHor[:iStartIteration]
                vVer = vVer[:iStartIteration]
                height = height[:iStartIteration]
                gamma = gamma[:iStartIteration]
                time = time[:iStartIteration]
                vLim = vLim[:iStartIteration]

                alpha.extend(alphaFlare)
                vHor.extend(vHorFlare)
                vVer.extend(vVerFlare)
                height.extend(heightFlare)
                gamma.extend(gammaFlare)
                time.extend(timeFlare)
                vLim.extend(vLimFlare)

                break

            # If this position is reached then a new starting height is set, but
            # the bias maps still need resetting
            biasBaseFlareGamma = defaultBiasBaseGamma
            biasBaseFlareVInf = defaultBiasBaseVInf
            biasBaseFlareGammaDot = defaultBiasBaseGammaDot
            biasFlareGamma.reset(biasBaseFlareGamma)
            biasFlareVInf.reset(biasBaseFlareVInf)
            biasFlareGammaDot.reset(biasBaseFlareGammaDot)

    # Rerun simulation to see if averaging out the initial values yields
    # approximately the same results
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
            alphaMin = int(max(iIteration - int(5.5 / dt), 0))
            alphaMax = int(min(iIteration + int(5.5 / dt), len(alpha)))

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
        # TODO: SAVE INTERPOLATED RESULTS AS WELL
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
        file.save('dive_' + str(heightUpper) + 'to' + str(heightTarget) +
                  '_' + str(vHorInitial) + 'to' + str(vHorTarget) +
                  '_' + str(vVerInitial) + 'to' + str(vVerTarget) + '.dat')

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
    #OptimizeDive(54999, 30000, -20, 10, 0, 0, 700*8.8, 35.0, 0, 0.01, lookupCl, lookupCd)
    #OptimizeDive(55000, 38000, -20, -10, 0, 0, 700*8.8, 35.0, 0, 0.25, lookupCl, lookupCd)
    #OptimizeDive(55000, 46000, -20, -10, 0, 0, 700*8.8, 35.0, 0, 0.25, lookupCl, lookupCd)
    OptimizeDive(55000, 35000, -20, 0, 0, 0, 700*8.8, 35.0, 20, -20, 0.15, lookupCl, lookupCd, storeResults=True)
    #OptimizeDive(38000, 30000, 50, 0, 0, 0, 700*8.8, 35.0, 0, 0.25, lookupCl, lookupCd)
    #OptimizeDive(46000, 30000, 35, -10, 0, 0, 700*8.8, 35.0, 0, 0.25, lookupCl, lookupCd)

def __TestPlotDive__():
    PlotDive("dive55000-30000.-20-0.dat")

__TestOptimizeDive__()
#__TestPlotDive__()
