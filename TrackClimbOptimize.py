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
import TimeEstimator
import TrackSettings

def OptimizeClimb(heightLower, heightUpper, heightQuit, vHorInitial, vVerInitial,
                  longitude, latitude, W, S, vHorTarget, vVerTarget, PRequired,
                  speedOfSoundRatio, inclination, dt, lookupCl, lookupCd,
                  severity=0.0, lookupBoundLowerVInf=None,
                  lookupBoundUpperVInf=None, plotResults=True, storeResults=True):
    # Construct lookup tables and interpolators
    atmosphere = Atmosphere.Atmosphere()
    lookupdCldAlpha = lookupCl.getDerivative()
    lookupdCddAlpha = lookupCd.getDerivative()
    lookupReverseCl = TrackLookup.LookupSegmented1D(lookupCl.getPoints()[1],
        lookupCl.getPoints()[0])

    if (not isinstance(PRequired, TrackLookup.Lookup1D)) and \
            (not isinstance(PRequired, TrackLookup.LookupSegmented1D)):
        # Convert the (hopefully scalar, otherwise it is error time) required
        # power into a lookup table
        offset = heightLower - heightQuit
        PRequired = TrackLookup.Lookup1D([heightQuit, heightUpper + offset],
            [PRequired, PRequired])


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
    biasWidth = 5000 # m

    alphaDotLimit = 0.5 # degree/s
    gammaDotLimit = 0.2 / 180.0 * np.pi # rad/s
    gammaLimit = np.pi / 2.0 # rad

    climboutFailHeight = heightUpper + (heightLower - heightQuit) * 0.1
    climboutGammaValid = 0.5 / 180.0 * np.pi # one-side of a two-sided range of gamma validity after climbout
    climboutHeightValid = 125 # one-side of a two-sided range of height validity after climbout

    updateCount = 35 # number of iterations before printing an update statement
    averageTime = 2.5 # number of seconds to average from the results for the resimulation
    disregardDotsTime = 5.0 # number of seconds to disregard alphaDot and gammaDot limits

    # Set initial values
    print(TrackCommon.StringHeader("Optimizing Climbing", 60))

    initialRho = TrackCommon.AdjustSeverity(atmosphere.density(heightLower, latitude, longitude), severity)
    initialZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(heightLower, latitude, longitude), severity)
    initialVInf = np.sqrt(np.power(initialZonal + vHorInitial, 2.0) + np.power(vVerInitial, 2.0))
    initialGamma = np.arctan2(vVerInitial, initialZonal + vHorInitial)

    initialAlpha = 0.0

    if abs(initialGamma) < 0.1 / 180.0 * np.pi:
        initialAlpha = TrackAngleOfAttack.AngleOfAttackSteady(W, S,
            0.5 * initialRho * initialVInf**2.0, lookupReverseCl)
    else:
        print(' > initial PReq:', PRequired(heightLower))
        initialAlpha = TrackAngleOfAttack.AngleOfAttackPowered(W, S,
            0.5 * initialRho * initialVInf**2.0, PRequired(heightLower) / initialVInf,
            initialGamma, inclination, lookupCl, lookupdCldAlpha)

    if initialAlpha[1] == False:
        raise ValueError("Initial angle of attack is invalid")

    initialAlpha = initialAlpha[0]
    print(' > initial alpha:', initialAlpha, 'degrees')
    print(' > initial gamma:', initialGamma * 180.0 / np.pi, 'degrees')
    print(' > initial PReq :', PRequired(heightLower) / 1e3, 'kW')
    #return

    alpha = [initialAlpha]
    vHor = [vHorInitial]
    vVer = [vVerInitial]
    height = [heightLower]
    gamma = [initialGamma]
    power = [PRequired(heightLower)]
    time = [0.0]
    vLim = [0.0]

    # Set bias maps and their associated variables
    biasBaseDefaultGamma = 0.975
    biasBaseDefaultVInf = 0.975
    biasBaseDefaultGammaDot = 0.975
    biasBaseDefaultBoundVInf = 0.975
    biasBaseDefaultVPositive = 0.975

    biasBaseGamma = biasBaseDefaultGamma
    biasBaseVInf = biasBaseDefaultVInf
    biasBaseGammaDot = biasBaseDefaultGammaDot
    biasBaseBoundLowerVInf = biasBaseDefaultBoundVInf
    biasBaseBoundUpperVInf = biasBaseDefaultBoundVInf
    biasBaseVPositive = biasBaseDefaultVPositive

    heightBiasLower = heightQuit
    heightBiasUpper = heightUpper + (heightLower - heightQuit)

    biasGamma = TrackBiasMap.BiasMap("gamma", heightBiasLower, heightBiasUpper, 1024, biasBaseGamma)
    biasVInf = TrackBiasMap.BiasMap("vInf", heightBiasLower, heightBiasUpper, 1024, biasBaseVInf)
    biasGammaDot = TrackBiasMap.BiasMap("gammaDot", heightBiasLower, heightBiasUpper, 1024, biasBaseGammaDot)
    biasBoundLowerVInf = TrackBiasMap.BiasMap("vInfLower", heightBiasLower, heightBiasUpper, 1024, biasBaseBoundLowerVInf)
    biasBoundUpperVInf = TrackBiasMap.BiasMap("vInfUpper", heightBiasLower, heightBiasUpper, 1024, biasBaseBoundUpperVInf)
    biasVPositive = TrackBiasMap.BiasMap("vPositive", heightBiasLower, heightBiasUpper, 1024, biasBaseVPositive)

    # Start iterating
    solved = False
    failed = False

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
                power[-1], W, S, inclination, alphaNew, dt, lookupCl, lookupCd,
                atmosphere, severity, tol=1e-8, relax=0.8)

            totalTime = dt * (iIteration + 1)

            # Filter the valid solutions from the invalid ones. Keep track of
            # the number of offenders to know how to alter the bias maps
            iValid = []
            vZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(hNew, latitude, longitude), severity)
            vInf = np.sqrt(np.power(vHorNew + vZonal, 2.0) + np.power(vVerNew, 2.0))
            vLimit = atmosphere.speedOfSound(hNew, latitude, longitude) * speedOfSoundRatio

            gammaOffenders = 0
            vInfOffenders = 0
            gammaDotOffenders = 0
            boundLowerVInfOffenders = 0
            boundUpperVInfOffenders = 0
            vPositiveOffenders = 0

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

                if iIteration >= int(disregardDotsTime / dt):
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

                if vHorNew[i] + vZonal[i] < 0:
                    vPositiveOffenders += 1
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

#                toShow = int(len(alphaNew) / 2)
#                print(TrackCommon.StringPad(" > gamma          = ", gammaNew[toShow] * 180.0 / np.pi, 3, 10) + " deg")
#                print(TrackCommon.StringPad(" > vInf           = ", vInf[toShow], 3, 10) + " m/s")
#                print(TrackCommon.StringPad(" > vLimit         = ", vLimit[toShow], 3, 10) + " m/s")
#                print(TrackCommon.StringPad(" > vZonal         = ", vZonal[toShow], 3, 10) + " m/s")
#                print(TrackCommon.StringPad(" > vHor           = ", vHorNew[toShow], 3, 10) + " m/s")
#                print(TrackCommon.StringPad(" > vVer           = ", vVerNew[toShow], 3, 10) + " m/s")
#                print(TrackCommon.StringPad(" > gammaDot       = ", abs((gammaNew[toShow] - gammaOld) * 180.0 / np.pi / dt), 5, 10) + " deg/s")
#                print(TrackCommon.StringPad(" > gammaDot limit = ", gammaDotLimit * 180.0 / np.pi, 5, 10) + " deg/s")

#                if lookupBoundLowerVInf != None:
#                    print(TrackCommon.StringPad(" > vHor lower     = ", lookupBoundLowerVInf(hNew[toShow]), 3, 10) + " m/s")
#
#                if lookupBoundUpperVInf != None:
#                    print(TrackCommon.StringPad(" > vHor upper     = ", lookupBoundUpperVInf(hNew[toShow]), 3, 10) + " m/s")

                # Determine how to adjust the bias maps
                listOffenders = [gammaOffenders, vInfOffenders, gammaDotOffenders,
                    boundLowerVInfOffenders, boundUpperVInfOffenders]
                iWorstOffender = np.argmax(listOffenders)

                if listOffenders[iWorstOffender] == 0:
                    print('\n * No offenders, adjusting base biases:')

                    biasBaseGamma, biasBaseVInf, biasBaseGammaDot, biasBaseBoundLowerVInf, \
                        biasBaseBoundUpperVInf, biasBaseVPositive = \
                        TrackCommon.AdjustBiasMapCommonly([biasGamma, biasVInf,
                        biasGammaDot, biasBoundLowerVInf, biasBoundUpperVInf,
                        biasVPositive], biasStep, ['gamma', 'vInf', 'gammaDot',
                        'vLower', 'vUpper', 'vPositive'])

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
                    elif iWorstOffender == 5:
                        biasBaseVPositive = TrackCommon.AdjustBiasMapIndividually(
                            biasVPositive, biasStep, height[-1], biasWidth, 'vPositive')
                    else:
                        raise RuntimeError("Unrecognized offender index for bias map")

                    print('')

                if biasBaseGamma < biasLimit or biasBaseVInf < biasLimit or \
                    biasBaseGammaDot < biasLimit or biasBaseBoundLowerVInf < biasLimit or \
                    biasBaseBoundUpperVInf < biasLimit or biasBaseVPositive < biasLimit:
                    print("Failed to find a solution")
                    failed = True
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
                #     a constant deltaTime, is the same as maximizing deltaHeight)
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
                        #    TrackCommon.StringPad(' hor: ', vHorNew[i], 1, 6) + " m/s" +
                        #    TrackCommon.StringPad(', lim: ', curBoundUpperVInf, 1, 6) + " m/s" +
                        #    TrackCommon.StringPad(', div: ', divFactor, 1, 6) + " m/s" +
                        #    TrackCommon.StringPad(', metric: ', ((vHorNew[i] - curBoundUpperVInf) / divFactor - curBiasBoundUpperVInf)**2.0, 3, 6))

                        metric += ((vHorNew[i] + curBoundDelta - curBoundUpperVInf) / curBoundDelta)**2.0

                # - influence of going in the direction of the wind
                curBiasVPositive = 1.0 - biasVPositive(hNew[i])
                totalVHor = vHorNew[i] + vZonal[i]
                if totalVHor < curBiasVPositive * vLimit[i]:
                    metric += ((curBiasVPositive * vLimit[i] - totalVHor) /
                        (curBiasVPositive * vLimit[i]))**2.0

                # TODO: TEST DISTANCE FROM ALPHA CL/CD MAX
                #metric += ((alphaNew[i] - alphaMaxClCd) / ClPoints[0][-1] *
                #           vInf[i] / vLimit[i])**2.0

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
            power.append(PRequired(hNew[iSolution]))
            time.append(totalTime)

            gammaOld = gammaNew[iSolution]

#            print(round(totalTime, 1), 's, h:', round(hNew[iSolution], 1), 'm, alpha:',
#                  round(alphaNew[iSolution], 1), 'deg, gamma:',
#                  round(gammaNew[iSolution] * 180.0 / np.pi, 1), 'deg, hor:',
#                  round(vHorNew[iSolution], 2), 'm/s, ver:',
#                  round(vVerNew[iSolution], 2), 'm/s')
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

    # Start the climbout phase. Climbout starts halfway in the climb and is
    # performed using a kind of bisection algorithm
    print(TrackCommon.StringHeader("Optimizing Climbout", 60))

    # Calculate final flight properties
    vZonalFinal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(heightUpper, latitude, longitude), severity)
    vInfFinal = np.sqrt(np.power(vHorTarget + vZonalFinal, 2.0) + np.power(vVerTarget, 2.0))
    gammaFinal = np.arctan2(vVerTarget, vZonalFinal + vHorTarget)

    # Setup heights for bisectional climbout iterations
    iterativeUpperHeight = heightUpper
    iterativeLowerHeight = heightLower
    iterativeCenterHeight = (heightUpper + heightLower) / 2.0

    # Setup special bias maps for the climbout
    biasBaseClimboutGamma = biasBaseDefaultGamma
    biasBaseClimboutVInf = biasBaseDefaultVInf
    biasBaseClimboutGammaDot = biasBaseDefaultGammaDot
    biasBaseClimboutVPositive = biasBaseDefaultVPositive

    biasClimboutHeightLower = iterativeLowerHeight - 0.1 * (iterativeUpperHeight - iterativeLowerHeight)
    biasClimboutHeightUpper = iterativeUpperHeight + 0.1 * (iterativeUpperHeight - iterativeLowerHeight)

    biasClimboutGamma = TrackBiasMap.BiasMap("gamma", biasClimboutHeightLower, biasClimboutHeightUpper, 1024, biasBaseClimboutGamma)
    biasClimboutVInf = TrackBiasMap.BiasMap("vInf", biasClimboutHeightLower, biasClimboutHeightUpper, 1024, biasBaseClimboutVInf)
    biasClimboutGammaDot = TrackBiasMap.BiasMap("gammaDot", biasClimboutHeightLower, biasClimboutHeightUpper, 1024, biasBaseClimboutGammaDot)
    biasClimboutVPositive = TrackBiasMap.BiasMap("vPositive", biasClimboutHeightLower, biasClimboutHeightUpper, 1024, biasBaseClimboutVPositive)

    print(TrackCommon.StringPad(" * Final height = ", heightUpper, 1, 10) + ' km')
    print(TrackCommon.StringPad(" * Final gamma  = ", gammaFinal * 180.0 / np.pi, 3, 10) + " deg")
    print(TrackCommon.StringPad(" * Target vInf  = ", vInfFinal, 2, 10) + " m/s")

    while iterativeCenterHeight < climboutFailHeight and (not failed):
        # Find which iteration step this height corresponds to
        iStartIteration = 0

        for i in range(0, len(height)):
            if height[i] > iterativeCenterHeight:
                if i == 0:
                    raise ValueError("Initial height is larger than center height. " +
                        "This should be impossible!")

                iStartIteration = i - 1
                break

        print(TrackCommon.StringPad("\n * Begin climb = ", height[iStartIteration], 1, 10) + " m\n")

        # Setup initial values for climbout
        alphaClimbout = [alpha[iStartIteration]]
        vHorClimbout = [vHor[iStartIteration]]
        vVerClimbout = [vVer[iStartIteration]]
        heightClimbout = [height[iStartIteration]]
        gammaClimbout = [gamma[iStartIteration]]
        timeClimbout = [time[iStartIteration]]
        vLimClimbout = [vLim[iStartIteration]]
        powerClimbout = [power[iStartIteration]]

        gammaOld = gammaClimbout[0]
        iIteration = 0
        solved = True

        while heightClimbout[-1] < heightUpper + climboutHeightValid and \
                abs(gammaClimbout[-1] - gammaFinal) > climboutGammaValid:
            # Determine new valid range of angles of attack
            alphaNew = np.linspace(max(alphaLimits[0], alphaClimbout[-1] - dt * alphaDotLimit),
                                   min(alphaLimits[1], alphaClimbout[-1] + dt * alphaDotLimit), numAlpha)

            # Determine new flight state
            vHorNew, vVerNew, gammaNew, hNew = TrackClimb.Step(heightClimbout[-1],
                alphaClimbout[-1], gammaClimbout[-1], vHorClimbout[-1],
                vVerClimbout[-1], longitude, latitude, powerClimbout[-1], W, S,
                inclination, alphaNew, dt, lookupCl, lookupCd, atmosphere,
                severity, tol=1e-8, relax=0.8)

            totalTime = timeClimbout[0] + dt * (iIteration + 1)

            # Filter the valid solutions from the invalid ones
            iValid = []
            vZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(hNew, latitude, longitude), severity)
            vInf = np.sqrt(np.power(vHorNew + vZonal, 2.0) + np.power(vVerNew, 2.0))
            vLimit = atmosphere.speedOfSound(hNew, latitude, longitude) * speedOfSoundRatio
            gammaDot = (gammaNew - gammaOld) / dt

            gammaOffenders = 0
            vInfOffenders = 0
            gammaDotOffenders = 0
            vPositiveOffenders = 0

            for i in range(0, len(alphaNew)):
                # Determine if the solution is valid
                isOffender = False

                if gammaNew[i] < -gammaLimit or gammaNew[i] > gammaLimit:
                    gammaOffenders += 1
                    isOffender = True

                if vInf[i] > vLimit[i]:
                    vInfOffenders += 1
                    isOffender = True

                if iIteration >= int(disregardDotsTime / dt):
                    if abs((gammaNew[i] - gammaOld) / dt) > gammaDotLimit:
                        gammaDotOffenders += 1
                        isOffender = True

                if vHorNew[i] + vZonal[i] < 0:
                    vPositiveOffenders += 1
                    isOffender = True

                if isOffender:
                    continue

                # The current values provide a valid solution
                iValid.append(i)

            if len(iValid) == 0:
                # No valid solutions were found
                print('\n * Did not find a solution at:\n > ' +
                    TrackCommon.StringPad('t = ', totalTime, 3, 10) + ' s\n > ' +
                    TrackCommon.StringPad('h = ', hNew[-1], 3, 10) + ' m\n')

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
                listOffenders = [gammaOffenders, vInfOffenders,
                    gammaDotOffenders, vPositiveOffenders]
                iWorstOffender = np.argmax(listOffenders)

                if listOffenders[iWorstOffender] == 0:
                    print('\n * No offenders, adjusting base biases:')
                    biasBaseClimboutGamma, biasBaseClimboutVInf,
                    biasBaseClimboutGammaDot, biasBaseClimboutVPositive = \
                        TrackCommon.AdjustBiasMapCommonly([biasBaseClimboutGamma,
                        biasBaseClimboutVInf, biasBaseClimboutGammaDot,
                        biasBaseClimboutVPositive], biasStep, ['gamma', 'vInf',
                        'gammaDot', 'vPositive'])
                    print('')
                else:
                    # Figure out how to adjust the bias map
                    print(TrackCommon.StringPad('\n * Adjusting bias, average gamma: ',
                        np.average(abs(gammaNew)) * 180.0 / np.pi, 2, 10))

                    if iWorstOffender == 0:
                        biasBaseClimboutGamma = TrackCommon.AdjustBiasMapIndividually(
                            biasClimboutGamma, biasStep, heightClimbout[-1], biasWidth, 'gamma')
                    elif iWorstOffender == 1:
                        biasBaseClimboutVInf = TrackCommon.AdjustBiasMapIndividually(
                            biasClimboutVInf, biasStep, heightClimbout[-1], biasWidth, 'vInf')
                    elif iWorstOffender == 2:
                        biasBaseClimboutGammaDot = TrackCommon.AdjustBiasMapIndividually(
                            biasClimboutGammaDot, biasStep, heightClimbout[-1], biasWidth, 'gammaDot')
                    elif iWorstOffender == 3:
                        biasBaseClimboutVPositive = TrackCommon.AdjustBiasMapIndividually(
                            biasClimboutVPositive, biasStep, heightClimbout[-1], biasWidth, 'vPositive')
                    else:
                        raise RuntimeError("Unrecognized climbout offender index for bias map")

                    print('')

                if biasBaseClimboutGamma < biasLimit or biasBaseClimboutVInf < biasLimit or \
                        biasBaseClimboutGammaDot < biasLimit or biasBaseClimboutVPositive < biasLimit:
                    print("Failed to find a climbout solution")
                    return

                # Restart with a new bias
                solved = False
                break

            # Seek the solution that minimizes the metric
            bestMetric = 1e19
            iSolution = 0

            for i in iValid:
                # Base metric contribution is to close in on the desired speed
                # and flight path angle. The distance to the final intended
                # speed is scaled linearly
                metric = ((vInf[i] - vInfFinal) / vLimit[i])**2.0
                metric += ((gammaNew[i] - gammaFinal) / gammaLimit)**2.0

                # Modifying contributions to steer away from the undesired
                # regions
                curBiasGammaDot = biasClimboutGammaDot(hNew[i])
                if abs(gammaDot[i]) > curBiasGammaDot * gammaDotLimit:
                    metric += ((gammaDot[i] - gammaDotLimit * curBiasGammaDot) /
                        (gammaDotLimit * (1.0 - curBiasGammaDot)))**2.0

                curBiasVInf = biasClimboutVInf(hNew[i])
                if vInf[i] > vLimit[i] * curBiasVInf:
                    metric += ((vInf[i] - vLimit[i] * curBiasVInf) /
                        (vLimit[i] * (1.0 - curBiasVInf)))**2.0

                curBiasGamma = biasClimboutGamma(hNew[i])
                if abs(gammaNew[i]) > curBiasGamma * gammaLimit:
                    metric += ((abs(gammaNew[i]) - curBiasGamma * gammaLimit) /
                        (gammaLimit * (1.0 - curBiasGamma)))**2.0

                curBiasVPositive = 1.0 - biasVPositive(hNew[i])
                totalVHor = vZonal[i] + vHorNew[i]
                if totalVHor < curBiasVPositive * vLimit[i]:
                    metric += ((curBiasVPositive * vLimit[i] - totalVHor) /
                        (curBiasVPositive * vLimit[i]))**2.0

                if metric < bestMetric:
                    bestMetric = metric
                    iSolution = i

            # Append the new optimum-metric solution to the solution arrays
            alphaClimbout.append(alphaNew[iSolution])
            vHorClimbout.append(vHorNew[iSolution])
            vVerClimbout.append(vVerNew[iSolution])
            heightClimbout.append(hNew[iSolution])
            gammaClimbout.append(gammaNew[iSolution])
            timeClimbout.append(totalTime)
            vLimClimbout.append(vLimit[iSolution])
            powerClimbout.append(PRequired(hNew[iSolution]))

            # Check if the current solution adheres to the stipulated requirements
            if (abs(heightClimbout[-1] - heightUpper) < climboutHeightValid and \
                    abs(gammaClimbout[-1] - gammaFinal) < climboutGammaValid):
                # Found a solution
                print('\n * Solution found at:\n > ' +
                    TrackCommon.StringPad('t     = ', totalTime, 3, 10) + ' s\n > ' +
                    TrackCommon.StringPad('h     = ', heightClimbout[-1], 1, 10) + ' m\n > ' +
                    TrackCommon.StringPad('gamma = ', gammaClimbout[-1] * 180.0 / np.pi, 2, 10) + " deg\n")
                break

            gammaOld = gammaNew[iSolution]

            if iIteration % updateCount == 0:
                print(TrackCommon.StringPad("Solved at t = ", totalTime, 3, 8) +
                      TrackCommon.StringPad(" s, h = ", hNew[iSolution], 0, 7) +
                      TrackCommon.StringPad(" m, Vver = ", vVerNew[iSolution], 2, 6) +
                      TrackCommon.StringPad(" m/s, gamma = ", gammaNew[iSolution] * 180.0 / np.pi, 3, 8) +
                      TrackCommon.StringPad(" deg, alpha = ", alphaNew[iSolution], 3, 8) + " deg")

            iIteration += 1

        if solved == True:
            if abs(gammaClimbout[-1] - gammaFinal) > climboutGammaValid:
                # Climbout angle differs by too much
                if gammaClimbout[-1] < gammaFinal:
                    # Flight path angle is less steep. This will probably never
                    # happen (height should be reached first)
                    raise RuntimeError("Flight path angle was unexpectedly less " +
                        "steep when the target height was reached")
                else:
                    # Flight path angle is too steep, need to initiate climbout
                    # earlier
                    iterativeUpperHeight = iterativeCenterHeight
                    iterativeCenterHeight = (iterativeLowerHeight + iterativeCenterHeight) / 2.0
                    print('\n * Flight path angle too steep, initiating climbout at',
                        round(iterativeCenterHeight / 1e3, 3), 'km\n')
            elif abs(heightClimbout[-1] - heightUpper) > climboutHeightValid:
                # Intended flight path angle is reached but the height not yet.
                # This means the climbout must be performed further up
                iterativeLowerHeight = iterativeCenterHeight
                iterativeCenterHeight = (iterativeLowerHeight + iterativeUpperHeight) / 2.0
                print('\n * Height at which climbout ended is too high, initiating climbout at',
                    round(iterativeCenterHeight / 1e3, 3), 'km\n')
            else:
                # Good solution, strip original solution and append the climbout
                alpha = alpha[:iStartIteration]
                vHor = vHor[:iStartIteration]
                vVer = vVer[:iStartIteration]
                height = height[:iStartIteration]
                gamma = gamma[:iStartIteration]
                time = time[:iStartIteration]
                vLim = vLim[:iStartIteration]
                power = power[:iStartIteration]

                alpha.extend(alphaClimbout)
                vHor.extend(vHorClimbout)
                vVer.extend(vVerClimbout)
                height.extend(heightClimbout)
                gamma.extend(gammaClimbout)
                time.extend(timeClimbout)
                vLim.extend(vLimClimbout)
                power.extend(powerClimbout)

                break

            # If this position is reached then a new starting height is set, but
            # the bias maps need to be reset
            biasBaseClimboutGamma = biasBaseDefaultGamma
            biasBaseClimboutVInf = biasBaseDefaultVInf
            biasBaseClimboutGammaDot = biasBaseDefaultGammaDot
            biasClimboutGamma.reset(biasBaseClimboutGamma)
            biasClimboutVInf.reset(biasBaseClimboutVInf)
            biasClimboutGammaDot.reset(biasBaseClimboutGammaDot)

    # Rerun the simulation with values that are averaged over a prolonged period
    # of time to see if the results match
    alphaFinal = [alpha[0]]
    vHorFinal = [vHor[0]]
    vVerFinal = [vVer[0]]
    heightFinal = [height[0]]
    gammaFinal = [gamma[0]]
    timeFinal = [time[0]]
    vLimFinal = [vLim[0]]
    powerFinal = [power[0]]

    numSteps = int(averageTime / dt)

    for iIteration in range(1, len(alpha)):
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
            latitude, powerFinal[-1], W, S, inclination, np.asarray([alphaFinal[-1]]),
            dt, lookupCl, lookupCd, atmosphere, severity, tol=1e-8, relax=0.8)

        # Store results
        alphaFinal.append(alphaAverage)
        vHorFinal.append(vHorNew[0])
        vVerFinal.append(vVerNew[0])
        heightFinal.append(hNew[0])
        gammaFinal.append(gammaNew[0])
        timeFinal.append(time[iIteration])
        powerFinal.append(PRequired(hNew[0]))

    # prepare Vinf
    vInf = np.zeros([len(vVer)])

    for i in range(0, len(vVer)):
        vZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(height[i], latitude, longitude), severity)
        vInf[i] = np.sqrt(np.power(vZonal + vHor[i], 2.0) + np.power(vVer[i], 2.0))

    vInfFinal = np.zeros([len(vVerFinal)])

    for i in range(0, len(vVerFinal)):
        vZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(heightFinal[i], latitude, longitude), severity)
        vInfFinal[i] = np.sqrt(np.power(vZonal + vHorFinal[i], 2.0) + np.power(vVerFinal[i], 2.0))

    # Plot the results
    if plotResults == True:
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

        axVinf.plot(time, vInf, 'g')
        axVinf.plot(timeFinal, vInfFinal, 'r')
        axVinf.plot(time, vLim, 'g--')
        axVinf.grid(True)

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

        # Plot the required power
        fig = plt.figure()
        axPower = fig.add_subplot(111)
        axPower.plot(time, np.asarray(power) / 1e3, 'g')
        axPower.plot(timeFinal, np.asarray(powerFinal) / 1e3, 'r')
        axPower.set_xlabel('time [s]')
        axPower.set_ylabel('power [kW]')

    if storeResults:
        # Create storage object
        file = TrackStorage.DataStorage()

        # Store all directly simulated variables
        file.addVariable('time', time)
        file.addVariable('alpha', alpha)
        file.addVariable('gamma', gamma)
        file.addVariable('height', height)
        file.addVariable('vVer', vVer)
        file.addVariable('vHor', vHor)
        file.addVariable('vInf', vInf)
        file.addVariable('power', power)
        file.addVariable('biasGamma', biasGamma.getMap(), [biasGamma.getAxis()])
        file.addVariable('biasVInf', biasVInf.getMap(), [biasVInf.getAxis()])
        file.addVariable('biasGammaDot', biasGammaDot.getMap(), [biasGammaDot.getMap()])
        file.addVariable('biasBoundLowerVInf', biasBoundLowerVInf.getMap(), [biasBoundLowerVInf.getAxis()])
        file.addVariable('biasBoundUpperVInf', biasBoundUpperVInf.getMap(), [biasBoundUpperVInf.getAxis()])
        file.addVariable('dt', dt)
        file.addVariable('vLim', vLim)
        file.addVariable('updateCount', updateCount)
        file.addVariable('averageTime', averageTime)

        # Store all variables post-simulated
        file.addVariable('timeFinal', timeFinal)
        file.addVariable('alphaFinal', alphaFinal)
        file.addVariable('gammaFinal', gammaFinal)
        file.addVariable('heightFinal', heightFinal)
        file.addVariable('vVerFinal', vVerFinal)
        file.addVariable('vHorFinal', vHorFinal)
        file.addVariable('vInfFinal', vInfFinal)
        file.addVariable('powerFinal', powerFinal)
        file.addVariable('vLimFinal', vLimFinal)

        # Save added variables to a file
        file.save('climb_' + str(heightLower) + 'to' + str(heightUpper) +
            '_' + str(round(vHorInitial, 1)) +
            '_' + str(round(vVerInitial, 1)) + '.dat')

    return timeFinal, heightFinal, vHorFinal, vVerFinal, vInfFinal, alphaFinal, gammaFinal

def PlotClimb(filename):
    # Load data from file
    file = TrackStorage.DataStorage()
    file.load(filename)

    # Retrieve original-simulation variables
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

    # Retrieve post-simulation variables
    timeFinal = file.getVariable("timeFinal").getValues()
    alphaFinal = file.getVariable("alphaFinal").getValues()
    gammaFinal = file.getVariable("gammaFinal").getValues()
    heightFinal = file.getVariable("heightFinal").getValues()
    vVerFinal = file.getVariable("vVerFinal").getValues()
    vHorFinal = file.getVariable("vHorFinal").getValues()
    powerFinal = file.getVariable("powerFinal").getValues()
    vLimFinal = file.getVariable("vLimFinal").getValues()

    fig = plt.figure()
    axSpeedDir = fig.add_subplot(221)
    axSpeed = axSpeedDir.twinx()
    axHeight = fig.add_subplot(222)
    axAnglesLeft = fig.add_subplot(223)
    axAnglesRight = axAnglesLeft.twinx()
    axAnglesDot = fig.add_subplot(224)

    # Plot all relevant speeds
    lSpeedVer, = axSpeed.plot(time, vVer, 'r', label=r'$V_{\mathrm{ver},i}$', alpha=0.3, linewidth=3)
    lSpeedHor, = axSpeed.plot(time, vHor, 'g', label=r'$V_{\mathrm{hor},i}$', alpha=0.3, linewidth=3)
    lSpeedInf, = axSpeed.plot(time, vInf, 'b', label=r'$V_{\mathrm{\infty}}$', alpha=0.3, linewidth=3)
    lSpeedLim, = axSpeed.plot(time, vLimit, 'k--', label=r'$V_{\mathrm{lim}}$')

    lSpeedFinalVer, = axSpeed.plot(timeFinal, vVerFinal, 'r', r'$V_{\mathrm{ver},f}$')
    lSpeedFinalHor, = axSpeed.plot(timeFinal, vHorFinal, 'g', r'$V_{\mathrm{hor},f}$')

    axSpeed.set_xlabel(r'$t\;[s]$')
    axSpeed.set_ylabel(r'$V\;[m/s]$')
    axSpeed.grid(True)
    axSpeed.legend()

    # Plot the height
    axHeight.plot(time, height / 1e3, 'r', label=r'$h_i$', alpha=0.3, linewidth=3)
    axHeight.plot(timeFinal, heightFinal / 1e3, 'r', label=r'$h_f$')

    axHeight.set_xlabel(r'$t\;[s]$')
    axHeight.set_ylabel(r'$h\;[km]$')
    axHeight.grid(True)

    # Plot the angles
    axAnglesLeft.plot(time, alpha, 'r', label=r'$\alpha_i$', alpha=0.3, linewidth=3)
    axAnglesLeft.plot(timeFinal, alphaFinal, 'r', label=r'$\alpha_f$')

    for tick in axAnglesLeft.get_yticklabels():
        tick.set_color('r')

    axAnglesRight.plot(time, gamma, 'g', label=r'$\gamma$', alpha=0.3, linewidth=3)
    axAnglesRight.plot(timeFinal, gammaFinal, 'g', label=r'')
    for tick in axAnglesRight.get_yticklabels():
        tick.set_color('g')

    axAnglesLeft.set_xlabel(r'$t\;[s]$')
    axAnglesLeft.set_ylabel(r'$\alpha\;[\degree]$')
    axAnglesRight.set_ylabel(r'$\gamma\;[\degree]$')

    # Plot the delta angles
    alphaDot = (alpha[1:] - alpha[0:-1]) / dt
    gammaDot = (gamma[1:] - gamma[0:-1]) / dt * 180.0 / np.pi
    axAnglesDot.plot(time[0:-1], alphaDot, 'r', label=r'$\mathrm{d}\alpha/\mathrm{d}t$', alpha=0.3, linewidth=3)
    axAnglesDot.plot(time[0:-1], gammaDot, 'g', label=r'$\mathrm{d}\gamma/\mathrm{d}t$', alpha=0.3, linewidth=3)
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

def GenerateAscentMaps(axisHeight, axisDeltaV, W, S, inclination, lookupCl,
                       lookupCd, atm, qInfMin, qInfMax, alphaMin, alphaMax,
                       PReqMin, PReqMax, severity=0.0,
                       storeResults=True):
    # Create derivatives of lookup maps
    lookupdCldAlpha = lookupCl.getDerivative()
    lookupdCddAlpha = lookupCd.getDerivative()
    axisVVer = np.linspace(0.1, 10.0, 50)

    # Settings for minima and maximum for minimum power estimation
    minDeltaVClearance = 5.0    # minimum velocity clearance from dangerous zones
    deltaVRadius = 5.0          # radius to construct around optimum ascent line

    # Preallocate result maps
    vVerMap = np.zeros([len(axisHeight), len(axisDeltaV)])
    PReqMap = np.zeros(vVerMap.shape)
    alphaMap = np.zeros(vVerMap.shape)
    qInfMap = np.zeros(vVerMap.shape)
    occludeMap = np.zeros(vVerMap.shape)

    # Retrieve atmospheric data and set up time estimator
    vZonal = TrackCommon.AdjustSeverity(atm.velocityZonal(axisHeight, 0, 0), severity)
    density = TrackCommon.AdjustSeverity(atm.density(axisHeight, 0, 0), severity)

    timeEstimator = TimeEstimator.TimeEstimator(len(axisHeight))
    timeEstimator.startTiming()

    # Keep track of when the first valid ascent value and the last valid
    # ascent value appears
    iFirstValid = -1
    iLastValid = -1
    stopValid = False

    print(TrackCommon.StringHeader("Generating Ascent Maps", 60))

    for iHeight in range(0, len(axisHeight)):
        timeEstimator.startIteration(iHeight)

        hasValid = False

        for iDeltaV in range(0, len(axisDeltaV)):
            vTotal = vZonal[iHeight] + axisDeltaV[iDeltaV]

            vVerOverPBest = -1e20
            vVerBest = 0
            PReqBest = 0
            alphaBest = 0
            qInfBest = 0
            found = False

            for iVVer in range(0, len(axisVVer)):
                vInfSquared = np.power(vTotal, 2.0) + np.power(axisVVer[iVVer], 2.0)
                qInf = 0.5 * density[iHeight] * vInfSquared
                gamma = np.arctan2(axisVVer[iVVer], vTotal)

                alpha, T, valid = TrackAngleOfAttack.AngleOfAttackThrustClimbing(
                    W, S, qInf, gamma, inclination, lookupCl, lookupdCldAlpha,
                    lookupCd, lookupdCddAlpha)

                if valid:
                    PReq = T * np.sqrt(vInfSquared)

                    if PReq >= PReqMin and PReq <= PReqMax and \
                            alpha >= alphaMin and alpha <= alphaMax:
                        vVerOverP = axisVVer[iVVer] / PReq
                        found = True

                        if vVerOverP > vVerOverPBest:
                            vVerOverPBest = vVerOverP
                            vVerBest = axisVVer[iVVer]
                            PReqBest = PReq
                            alphaBest = alpha
                            qInfBest = qInf

                    if PReq >= PReqMin and PReq <= PReqMax and \
                            alpha >= alphaMin and alpha <= alphaMax and \
                            qInf >= qInfMin and qInf <= qInfMax:
                        hasValid = True

            # Store the best values in the maps
            vVerMap[iHeight, iDeltaV] = vVerBest
            PReqMap[iHeight, iDeltaV] = PReqBest
            alphaMap[iHeight, iDeltaV] = alphaBest
            qInfMap[iHeight, iDeltaV] = qInfBest

            if found == True:
                occludeMap[iHeight, iDeltaV] = 1.0

        if hasValid:
            if iFirstValid == -1:
                iFirstValid = iHeight

            if not stopValid:
                iLastValid = iHeight
        else:
            if iFirstValid != -1:
                stopValid = True

        timeEstimator.finishedIteration(iHeight)

        print('elapsed:', timeEstimator.getTotalElapsed(),
              ', remaining:', timeEstimator.getEstimatedRemaining())

    # Go through the created maps and figure out the best path for ascent
    pathDeltaV = np.zeros([iLastValid - iFirstValid + 1])
    pathMinDeltaV = np.zeros(pathDeltaV.shape)
    pathMaxDeltaV = np.zeros(pathDeltaV.shape)
    pathHeight = np.zeros(pathDeltaV.shape)
    pathPReq = np.zeros(pathDeltaV.shape)
    pathVVer = np.zeros(pathDeltaV.shape)

    for iHeight in range(iFirstValid, iLastValid + 1):
        # Find the optimum deltaV for the current height value
        localHeightIndex = iHeight - iFirstValid
        pathHeight[localHeightIndex] = axisHeight[iHeight]
        highestVVerOverPReq = -1e20
        iHighestDeltaV = -1
        iMinimumDeltaV = 0
        iMaximumDeltaV = 0
        inValid = False

        # Loop through all deltaV values. While searching for the optimum value
        # attempt to find the minimum and maximum operative ones as well
        for iDeltaV in reversed(range(0, len(axisDeltaV))):
            curAlpha = alphaMap[iHeight, iDeltaV]
            curPReq = PReqMap[iHeight, iDeltaV]
            curVVer = vVerMap[iHeight, iDeltaV]
            curQInf = qInfMap[iHeight, iDeltaV]

            if curAlpha >= alphaMin and curAlpha <= alphaMax and \
                    curPReq >= PReqMin and curPReq <= PReqMax and \
                    curQInf >= qInfMin and curQInf <= qInfMax:

                if not inValid:
                    # Found the first valid value
                    iMaximumDeltaV = iDeltaV
                    inValid = True

                # Check if this value is a new optimum
                vVerOverPReq = curVVer / curPReq

                if vVerOverPReq > highestVVerOverPReq:
                    iHighestDeltaV = iDeltaV
                    highestVVerOverPReq = vVerOverPReq
            elif inValid == True:
                # Exited the valid region
                iMinimumDeltaV = iDeltaV
                break

        # Iterate away until the desired clearance is achieved
        if axisDeltaV[iMaximumDeltaV] - axisDeltaV[iMinimumDeltaV] < 2.0 * minDeltaVClearance:
            # No index can be found where the optimum value is removed from the
            # valid boundaries by the given clearance, take the average
            iHighestDeltaV = int((iMinimumDeltaV + iMaximumDeltaV) / 2)
            pathDeltaV[localHeightIndex] = axisDeltaV[iHighestDeltaV]
            pathMinDeltaV[localHeightIndex] = axisDeltaV[iMinimumDeltaV]
            pathMaxDeltaV[localHeightIndex] = axisDeltaV[iMaximumDeltaV]
        else:
            while axisDeltaV[iHighestDeltaV] - axisDeltaV[iMinimumDeltaV] < minDeltaVClearance:
                iHighestDeltaV += 1

            while axisDeltaV[iMaximumDeltaV] - axisDeltaV[iHighestDeltaV] < minDeltaVClearance:
                iHighestDeltaV -= 1

            pathDeltaV[localHeightIndex] = axisDeltaV[iHighestDeltaV]
            pathMinDeltaV[localHeightIndex] = pathDeltaV[localHeightIndex] - deltaVRadius
            pathMaxDeltaV[localHeightIndex] = pathDeltaV[localHeightIndex] + deltaVRadius

        pathPReq[localHeightIndex] = PReqMap[iHeight, iHighestDeltaV]
        pathVVer[localHeightIndex] = vVerMap[iHeight, iHighestDeltaV]

    if storeResults:
        file = TrackStorage.DataStorage()
        file.addVariable('height', pathHeight)
        file.addVariable('minDeltaV', pathMinDeltaV)
        file.addVariable('avgDeltaV', pathDeltaV)
        file.addVariable('maxDeltaV', pathMaxDeltaV)
        file.addVariable('pReq', pathPReq)
        file.addVariable('vVer', pathVVer)

        file.save('optclimb_' + str(round(axisDeltaV[0], 3)) +
                  'to' + str(round(axisDeltaV[-1], 3)) +
                  '_' + str(round(severity, 2)) + '.dat')

    return {
        'axisHeight': axisHeight, 'axisDeltaV':axisDeltaV,
        'vVer': vVerMap, 'PReq':PReqMap, 'alpha': alphaMap,
        'qInf': qInfMap, 'occlude': occludeMap, 'pathHeight': pathHeight,
        'pathMinDeltaV': pathMinDeltaV, 'pathMaxDeltaV': pathMaxDeltaV,
        'pathVVer': pathVVer
    }

def PlotAscentMaps(axisHeight, axisDeltaV, vVerMap, PReqMap, alphaMap, qInfMap,
                   occludeMap, pathMinDeltaV, pathMaxDeltaV, pathHeight):
    # Plot the maps
    fig = plt.figure()
    bordersPRatioData, bordersPRatioLegend = TrackCommon.ImageAxes(0.0, 0.5, 0.5, 1.0)
    bordersPReqData, bordersPReqLegend = TrackCommon.ImageAxes(0.5, 1.0, 0.5, 1.0)
    bordersVVerData, bordersVVerLegend = TrackCommon.ImageAxes(0.0, 0.5, 0.0, 0.5)
    bordersAlphaData, bordersAlphaLegend = TrackCommon.ImageAxes(0.5, 1.0, 0.0, 0.5)
    occludeColor='k'

    absMinPReqRatio = 1e20
    absMaxPReqRatio = -1e20

    for i in range(0, len(axisHeight)):
        for j in range(0, len(axisDeltaV)):
            if abs(occludeMap[i, j] - 1.0) < 1e-1:
                curRatio = vVerMap[i, j] / PReqMap[i, j]

                if curRatio < absMinPReqRatio:
                    absMinPReqRatio = curRatio

                if curRatio > absMaxPReqRatio:
                    absMaxPReqRatio = curRatio

    # Plot the height to energy ratio data
    axPRatioData = fig.add_axes(bordersPRatioData)
    axPRatioLegend = fig.add_axes(bordersPRatioLegend)

    TrackCommon.PlotImage(fig, axPRatioData, axPRatioLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', vVerMap / (PReqMap / 1e6), r'$\Delta h / \Delta E\;[m/MJ]$',
        forceNormMin=absMinPReqRatio / 2 * 1e6, forceNormMax=absMaxPReqRatio * 1e6)

    axPRatioData.plot(pathMinDeltaV, pathHeight / 1e3, 'g--')
    axPRatioData.plot((pathMinDeltaV + pathMaxDeltaV) / 2.0, pathHeight / 1e3, 'g')
    axPRatioData.plot(pathMaxDeltaV, pathHeight / 1e3, 'g--')

    # Plot the power required data
    axPReqData = fig.add_axes(bordersPReqData)
    axPReqLegend = fig.add_axes(bordersPReqLegend)

    TrackCommon.PlotImage(fig, axPReqData, axPReqLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', qInfMap, r'$q_\infty\;[Pa]$')

    axPReqData.plot(pathMinDeltaV, pathHeight / 1e3, 'g--')
    axPReqData.plot((pathMinDeltaV + pathMaxDeltaV) / 2.0, pathHeight / 1e3, 'g')
    axPReqData.plot(pathMaxDeltaV, pathHeight / 1e3, 'g--')

    # Plot the vertical speed data
    axVVerData = fig.add_axes(bordersVVerData)
    axVVerLegend = fig.add_axes(bordersVVerLegend)

    TrackCommon.PlotImage(fig, axVVerData, axVVerLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', vVerMap, r'$V_\mathrm{ver}\;[m/s]$')

    axVVerData.plot(pathMinDeltaV, pathHeight / 1e3, 'g--')
    axVVerData.plot((pathMinDeltaV + pathMaxDeltaV) / 2.0, pathHeight / 1e3, 'g')
    axVVerData.plot(pathMaxDeltaV, pathHeight / 1e3, 'g--')

    # Plot the angle of attack data
    axAlphaData = fig.add_axes(bordersAlphaData)
    axAlphaLegend = fig.add_axes(bordersAlphaLegend)

    TrackCommon.PlotImage(fig, axAlphaData, axAlphaLegend, axisDeltaV, r'$\Delta V\;[m/s]$',
        axisHeight / 1e3, r'$h\;[km]$', alphaMap, r'$\alpha\;[\degree]$')

    axAlphaData.plot(pathMinDeltaV, pathHeight / 1e3, 'g--')
    axAlphaData.plot((pathMinDeltaV + pathMaxDeltaV) / 2.0, pathHeight / 1e3, 'g')
    axAlphaData.plot(pathMaxDeltaV, pathHeight / 1e3, 'g--')

    # Apply occlusion map to the various plots
    axPRatioData.contourf(axisDeltaV, axisHeight / 1e3, occludeMap,
                          [-0.5, 0.5], origin='lower', extend='min',
                          colors=[occludeColor, occludeColor])
    axPReqData.contourf(axisDeltaV, axisHeight / 1e3, occludeMap,
                        [-0.5, 0.5], origin='lower', extend='min',
                        colors=[occludeColor, occludeColor])
    axAlphaData.contourf(axisDeltaV, axisHeight / 1e3, occludeMap,
                         [-0.5, 0.5], origin='lower', extend='min',
                         colors=[occludeColor, occludeColor])

def __TestOptimizeClimb__():
    settings = TrackSettings.Settings()
    OptimizeClimb(38000, 62000, 30000, -10, 0, settings.longitude,
        settings.latitude, settings.W, settings.S, -10, 10, 40000, 0.7,
        settings.inclination, 0.25, settings.lookupCl, settings.lookupCd, 0.0,
        settings.lookupLower, settings.lookupUpper)

def __TestOptimizeBoundedClimb__(filename, lower, higher, finalVHor, finalVVer,
                                 initialVHor = None, initialVVer = None,
                                 severity = 0.0):
    settings = TrackSettings.Settings()

    file = TrackStorage.DataStorage()
    file.load(filename)

    axisHeight = file.getVariable('height').getValues()
    axisMin = file.getVariable('minDeltaV').getValues()
    axisMax = file.getVariable('maxDeltaV').getValues()
    vVer = file.getVariable('vVer').getValues()

    lookupVVer = TrackLookup.Lookup1D(axisHeight, vVer)

    if initialVHor == None:
        initialVHor = (settings.lowerBound(lower) + settings.upperBound(lower)) / 2.0

    if initialVVer == None:
        initialVVer = lookupVVer(lower)

    print('initial vHor =', round(initialVHor, 4), 'm/s')
    print('initial vVer =', round(initialVVer, 4), 'm/s')

    OptimizeClimb(lower, higher, 30000, initialVHor, initialVVer, 0, 0,
        settings.W, settings.S, finalVHor, finalVVer, 32e3,
        settings.speedOfSoundRatio, settings.inclination, 0.25,
        settings.lookupCl, settings.lookupCd, severity, settings.lowerBound,
        settings.upperBound)

def __TestAscentMap__(severity, vMin, vMax, qInfMin, qInfMax):
    settings = TrackSettings.Settings()
    atm = Atmosphere.Atmosphere()
    axisHeight = np.linspace(15, 80, 150) * 1e3
    axisDeltaV = np.linspace(vMin, vMax, 150)

    returned = GenerateAscentMaps(axisHeight, axisDeltaV,
        settings.W, settings.S, settings.inclination, settings.lookupCl,
        settings.lookupCd, atm, qInfMin, qInfMax, -8.0, 8.0, 0, 32e3, severity,
        storeResults=False)

    axisHeight = returned['axisHeight']
    axisDeltaV = returned['axisDeltaV']
    vVerMap = returned['vVer']
    pReqMap = returned['PReq']
    alphaMap = returned['alpha']
    qInfMap = returned['qInf']
    occludeMap = returned['occlude']
    pathHeight = returned['pathHeight']
    pathMinDeltaV = returned['pathMinDeltaV']
    pathMaxDeltaV = returned['pathMaxDeltaV']
    pathVVer = returned['pathVVer']

    PlotAscentMaps(axisHeight, axisDeltaV, vVerMap, pReqMap, alphaMap, qInfMap,
        occludeMap, pathMinDeltaV, pathMaxDeltaV, pathHeight)

#__TestOptimizeClimb__()
#PlotClimb('climb_35000to50000_50000_-10_0.dat')
#__TestOptimizeBoundedClimb__('./optclimb_-60.0to20.0_0.0.dat', 38000, 62000, 7.8, 0, -5.0, 0.0, 0.0)
#__TestAscentMap__(-1.6, -80, 5)
#__TestAscentMap__(0.0, -60, 20, 200, 1e10)
#__TestAscentMap__(-1.0, -60, 20, 200, 1e10)
