# -*- coding: utf-8 -*-
"""
Created on Sun May 29 12:43:41 2016

@author: MaxHenger
"""

import Atmosphere
import TrackLookup
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as scp_int

# Step is the function that performs a single iteration with regards to
# simulating the response while climibing under a given angle of attack at a
# constant required propellor power. The timestep should be relatively small as
# the low order integrations used can quickly cause instabilities in the
# algorithm.
# Input:
#   - hCur: The current height, a scalar in meters
#   - alphaCur: The current angle of attack, a scalar in degrees
#   - gammaCur: The current climbing angle, a scalar in radians
#   - vHorCur = The current horizontal velocity, a scalar in meters per second
#   - vVerCur = The current vertical velocity, a scalar in meters per second
#   - longitudeCur: The approximate solar longitude in degrees
#   - latitudeCur: The approximate latitude in degrees
#   - PRequired: The required propellor power, a scalar in Watt
#   - W: The weight of the aircraft, a scalar in Newton
#   - S: The surface area of the wing planform of the aircraft, a scalar in
#       meters squared
#   - inclination: The inclination of the propellors with respect to the
#       aircraft body, a scalar in radians
#   - alphaNew: The new angle of attack, an array in degrees
#   - dt: The timestep used in the simulation in seconds
#   - lookupCl: An instance of a TrackLookup class in which, for a given angle
#       of attack in degrees the lift coefficient can be obtained
#   - lookupCd: An instance of a TrackLookup class in which, for a given angle
#       of attack in degrees the drag coefficient can be obtained
#   - atmosphere: An instance of the Atmosphere class
#   - tol: The tolerance to use to determine when to stop iterating
def Step(hCur, alphaCur, gammaCur, vHorCur, vVerCur,
         longitudeCur, latitudeCur, PRequired, W, S, inclination,
         alphaNew, dt, lookupCl, lookupCd, atmosphere, tol=1e-8):
    # Calculate some often used variables
    # - general atmospheric variables
    rhoCur = atmosphere.density(hCur, latitudeCur, longitudeCur)[1]
    vWindZonalCur = atmosphere.velocityZonal(hCur, latitudeCur, longitudeCur)[1]
    gCur = 8.8 # TODO: Integrate with gravity model
    gNew = 8.8 # TODO: Integrate with gravity model
    vInfCurSquared = np.power(vHorCur + vWindZonalCur, 2.0) + np.power(vVerCur, 2.0)
    alphaNewRad = alphaNew / 180.0 * np.pi

    # - current variables related to generating lift, thrust and drag
    KA1 = PRequired * gCur / (W * np.sqrt(vInfCurSquared))
    KB1 = gCur * 0.5 * rhoCur * vInfCurSquared * S / W

    ClCur = lookupCl.find(alphaCur)
    ClCosCur = ClCur * np.cos(gammaCur)
    ClSinCur = ClCur * np.sin(gammaCur)

    CdCur = lookupCd.find(alphaCur)
    CdCosCur = CdCur * np.cos(gammaCur)
    CdSinCur = CdCur * np.sin(gammaCur)

    engineAngleCur = gammaCur + alphaCur / 180.0 * np.pi + inclination
    engineCosAngleCur = np.cos(engineAngleCur)
    engineSinAngleCur = np.sin(engineAngleCur)

    # - new variables related to generating lift, thrust and drag
    ClNew = np.zeros(alphaNew.shape)
    CdNew = np.zeros(alphaNew.shape)

    for iAlpha in range(0, len(alphaNew)):
        ClNew[iAlpha] = lookupCl.find(alphaNew[iAlpha])
        CdNew[iAlpha] = lookupCd.find(alphaNew[iAlpha])

    # Calculate the intital step's variables
    vHorNew = vHorCur + 0.5 * (
        KA1 * (engineCosAngleCur + np.cos(gammaCur + alphaNewRad + inclination)) -
        KB1 * ((ClSinCur + ClNew * np.sin(gammaCur)) + (CdCosCur + CdNew * np.cos(gammaCur)))
    ) * dt

    vVerNew = vVerCur + 0.5 * (
        KA1 * (engineSinAngleCur + np.sin(gammaCur + alphaNewRad + inclination)) +
        KB1 * ((ClCosCur + ClNew * np.cos(gammaCur)) - (CdSinCur + CdNew * np.sin(gammaCur))) -
        2 * gCur
    ) * dt

    # Start iterating until a stable solution is found
    for iIteration in range(0, 1000):
        # Calculate the new gamma and atmospheric properties
        hNew = hCur + 0.5 * (vVerCur + vVerNew) * dt
        vWindZonalNew = atmosphere.velocityZonal(hNew, latitudeCur, longitudeCur)[1]
        gammaNew = np.arctan2(vVerNew, vWindZonalNew + vHorNew)
        rhoNew = atmosphere.density(hNew, latitudeCur, longitudeCur)[1]
        vInfNewSquared = np.power(vHorNew + vWindZonalNew, 2.0) + np.power(vVerNew, 2.0)

        # Calculate new constants related to calculating acceleration
        KA2 = gNew * PRequired / (W * np.sqrt(vInfNewSquared))
        KB2 = gNew * 0.5 * rhoNew * vInfNewSquared * S / W

        # Calculate new horizontal and vertical speeds
        vHorNewGuess = vHorCur + 0.5 * (
            KA1 * engineCosAngleCur + KA2 * np.cos(gammaNew + alphaNewRad + inclination) -
            KB1 * (ClSinCur + CdCosCur) - KB2 * (ClNew * np.sin(gammaNew) + CdNew * np.cos(gammaNew))
        ) * dt

        vVerNewGuess = vVerCur + 0.5 * (
            KA1 * engineSinAngleCur + KA2 * np.sin(gammaNew + alphaNewRad + inclination) +
            KB1 * (ClCosCur - CdSinCur) + KB2 * (ClNew * np.cos(gammaNew) - CdNew * np.sin(gammaNew)) -
            (gCur + gNew)
        ) * dt

        # Check if the new estimations are all below the desired threshold
        vFactorVer = (vVerNewGuess - vVerNew) / vVerNew
        converged = True

        for iFactor in range(0, len(vFactorVer)):
            if vFactorVer[iFactor] > tol:
                # Vertical speed has not converged enough
                print('vFactorVer', iFactor, '=', vFactorVer[iFactor], 'too large ( iteration =', iIteration, ')')
                converged = False
                break

        if not converged:
            vHorNew = vHorNewGuess
            vVerNew = vVerNewGuess
            continue

        vFactorHor = (vHorNewGuess - vHorNew) / vHorNew

        for iFactor in range(0, len(vFactorHor)):
            if vFactorHor[iFactor] >tol:
                # Horizontal speed has not converged enough
                print('vFactorHor', iFactor, '=', vFactorHor[iFactor], 'too large ( iteration =', iIteration, ')')
                converged = False
                break

        if not converged:
            vHorNew = vHorNewGuess
            vVerNew = vVerNewGuess
            continue

        # If this position is reached then all values have converged to the
        # desired degree
        print('Solution converged after', iIteration + 1, 'iterations')
        return vHorNewGuess, vVerNewGuess, gammaNew

    raise RuntimeError("TrackClimb.Step did not converge")

def TestStep(h, alpha, gamma, alphaNew, inc, Vhor, Vver, dt, Pr, W, S, luts, atm):
    # Calculate often used variables
    rho = atm.density(h, 0, 0)
    vZonal = atm.velocityZonal(h, 0, 0)
    g1 = 8.8
    g2 = 8.8
    VinfSq = (Vhor + vZonal[1])**2.0 + Vver**2.0
    KA1 = Pr * g1 / (W * np.sqrt(VinfSq))
    KB1 = g1 * 0.5 * rho[1] * VinfSq * S / W

    Cl1 = luts[2].find(alpha)
    ClCos1 = Cl1 * np.cos(gamma)
    ClSin1 = Cl1 * np.sin(gamma)

    Cd1 = luts[4].find(alpha)
    CdCos1 = Cd1 * np.cos(gamma)
    CdSin1 = Cd1 * np.sin(gamma)

    Cl2 = np.zeros(alphaNew.shape)
    Cd2 = np.zeros(alphaNew.shape)

    for i in range(0, len(alphaNew)):
        Cl2[i] = luts[2].find(alphaNew[i])
        Cd2[i] = luts[4].find(alphaNew[i])

    alphaRad = alpha / 180.0 * np.pi
    alphaNewRad = alphaNew / 180.0 * np.pi

    # Calculate the initial Step
    VhorNew = Vhor + 0.5 * (
        KA1 * (np.cos(gamma + alphaRad + inc) + np.cos(gamma + alphaNewRad + inc)) -
        KB1 * ((ClSin1 + Cl2 * np.sin(gamma)) + (CdCos1 + Cd2 * np.cos(gamma)))
    ) * dt

    VverNew = Vver + 0.5 * (
        KA1 * (np.sin(gamma + alphaRad + inc) + np.sin(gamma + alphaNewRad + inc)) +
        KB1 * ((ClCos1 + Cl2 * np.cos(gamma)) - (CdSin1 + Cd2 * np.sin(gamma)))
        - 2 * g1
    ) * dt

    KA2 = np.zeros(alphaNew.shape)
    KB2 = np.zeros(alphaNew.shape)

    VhorAll = []
    VverAll = []
    gammaAll = []
    hAll = []
    VhorAll.append(VhorNew)
    VverAll.append(VverNew)

    for i in range(0, 10):
        hNew = h + 0.5 * (Vver + VverNew) * dt
        vZonalNew = atm.velocityZonal(hNew, 0, 0)[1]
        gammaNew = np.arctan2(VverNew, vZonalNew + VhorNew)
        rhoNew = atm.density(hNew, 0, 0)[1]
        VinfNewSq = np.power(VhorNew + vZonalNew, 2.0) + np.power(VverNew, 2.0)

        KA2 = g2 * Pr / (W * np.sqrt(VinfNewSq))
        KB2 = g2 * 0.5 * rhoNew * VinfNewSq * S / W

        VhorNew = Vhor + 0.5 * (
            KA1 * np.cos(gamma + alphaRad + inc) + KA2 * np.cos(gammaNew + alphaNewRad + inc) -
            KB1 * (ClSin1 + CdCos1) - KB2 * (Cl2 * np.sin(gammaNew) + Cd2 * np.cos(gammaNew))
        ) * dt

        VverNew = Vver + 0.5 * (
            KA1 * np.sin(gamma + alphaRad + inc) + KA2 * np.sin(gammaNew + alphaNewRad + inc) +
            KB1 * (ClCos1 - CdSin1) + KB2 * (Cl2 * np.cos(gammaNew) - Cd2 * np.sin(gammaNew)) -
            (g1 + g2)
        ) * dt

        VverAll.append(VverNew)
        VhorAll.append(VhorNew)
        gammaAll.append(gammaNew)
        hAll.append(hNew - h)

    VverAll = np.asarray(VverAll)
    VhorAll = np.asarray(VhorAll)
    gammaAll = np.asarray(gammaAll)
    hAll = np.asarray(hAll)

    fig = plt.figure()
    axVer = fig.add_subplot(221)
    axHor = fig.add_subplot(222)
    axGamma = fig.add_subplot(223)
    axH = fig.add_subplot(224)
    linesVer = []
    linesHor = []
    linesGamma = []
    linesH = []
    name = []
    nameGamma = []
    cmap = plt.get_cmap('jet')

    for i in range(0, VverAll.shape[0]):
        color = cmap(i / (VverAll.shape[0] - 1))
        lVer, = axVer.plot(alphaNew, VverAll[i], color=color)
        lHor, = axHor.plot(alphaNew, VhorAll[i], color=color)

        linesVer.append(lVer)
        linesHor.append(lHor)
        name.append('i = ' + str(i + 1))

    for i in range(0, gammaAll.shape[0]):
        color = cmap(i / (gammaAll.shape[0] - 1))
        lGamma, = axGamma.plot(alphaNew, gammaAll[i] * 180.0 / np.pi, color=color)
        lH, = axH.plot(alphaNew, hAll[i], color=color)
        linesGamma.append(lGamma)
        linesH.append(lH)
        nameGamma.append('i = ' + str(i + 1))

    axVer.legend(linesVer, name)
    axHor.legend(linesHor, name)
    axGamma.legend(linesGamma, nameGamma)
    axH.legend(linesH, nameGamma)
    axVer.set_title('vertical')
    axHor.set_title('horizontal')
    axGamma.set_title('gamma')
    axH.set_title('delta H')
    axVer.grid(True)
    axHor.grid(True)
    axGamma.grid(True)
    axH.grid(True)
    axVer.plot([alphaNew[0], alphaNew[-1]], [Vver, Vver], 'k--')
    axHor.plot([alphaNew[0], alphaNew[-1]], [Vhor, Vhor], 'k--')
    axGamma.plot([alphaNew[0], alphaNew[-1]], np.asarray([gamma, gamma]) * 180.0 / np.pi, 'k--')

#atm = Atmosphere.Atmosphere()
#luts = TrackDive.TestGenerateLiftOverDrag()
#TestStep(45000, 2.0, 2.0 / 180.0 * np.pi, np.linspace(-5, 5, 150), 0, 30, 4, 0.5,
#    80, 1, 0.001, luts, atm)
#Step(35000, 2.0, 2.0 / 180.0 * np.pi, 50, 5, 0, 0, 70, 1, 0.001, 0,
#    np.linspace(-5, 5, 50), 0.1, luts[2], luts[4], atm, tol=1e-15)
