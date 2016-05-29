# -*- coding: utf-8 -*-
"""
Created on Fri May 27 11:58:10 2016

@author: MaxHenger
"""

import Atmosphere
import TrackLookup
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as scp_int

# DiveStep is the function that performs a single iteration with regards to
# simulating the response while diving under a given angle of attack. The
# timestep should be relatively small as the low order integrations used can
# quickly cause instabilities in the algorithm.
# Input:
#   - hCur: The current height, a scalar in meters
#   - alphaCur: The current angle of attack, a scalar in degrees
#   - gammaCur: The current glide path angle, a scalar in radians
#   - vHorCur: The current horizontal velocity, a scalar in meters per second
#   - vVerCur: The current vertical velocity, a scalar in meters per second
#   - longitudeCur: The approximate solar longitude in degrees
#   - latitudeCur: The approximate latitude in degrees
#   - W: The weight of the aircraft, a scalar in Newton
#   - S: The surface area of the wing planform of the aircraft, a scalar in
#       meters squared
#   - alphaNew: The new angle of attack, an array in degrees
#   - dt: The timestep used in the simulation in seconds
#   - lookupCl: An instance of a TrackLookup class in which, for a given angle
#       of attack in degrees the lift coefficient can be obtained
#   - lookupCd: An instance of a TrackLookup class in which, for a given angle
#       of attack in degrees the drag coefficient can be obtained
#   - atmopshere: An instance of the Atmopshere class
#   - tol: The tolerance to use to determine when to stop iterating
#
#
# Output:
#   - vHorNew: New horizontal speeds, an array of the same length as alphaNew,
#       in meters per second
#   - vVerNew: New vertical speeds, an array of the same length as alphaNew,
#       in meters per second
#   - gammaNew: New glide path angles, an array of the same length as alphaNew,
#       in radians
def Step(hCur, alphaCur, gammaCur, vHorCur, vVerCur,
         longitudeCur, latitudeCur, W, S, alphaNew,
         dt, lookupCl, lookupCd, atmosphere, tol=1e-8):
    # Calculate some often used variables
    # - general atmospheric variables
    rhoCur = atmosphere.density(hCur, 0, 0)
    vWindZonalCur = atmosphere.velocityZonal(hCur, 0, 0)
    gCur = 8.8 # TODO: Update these with the gravity model
    gNew = 8.8 # TODO: Update these with the gravity model
    vInfCurSquared = (vHorCur + vWindZonalCur[1])**2.0 + vVerCur**2.0

    # - current variables related to generating lift and drag
    K1 = gCur * 0.5 * rhoCur[1] * vInfCurSquared * S / W

    ClCur = lookupCl.find(alphaCur)
    ClCosCur = ClCur * np.cos(gammaCur)
    ClSinCur = ClCur * np.sin(gammaCur)

    CdCur = lookupCd.find(alphaCur)
    CdCosCur = CdCur * np.cos(gammaCur)
    CdSinCur = CdCur * np.sin(gammaCur)

    # - new variables related to generating lift and drag
    ClNew = np.zeros(alphaNew.shape)
    CdNew = np.zeros(alphaNew.shape)

    for iAlpha in range(0, len(alphaNew)):
        ClNew[iAlpha] = lookupCl.find(alphaNew[iAlpha])
        CdNew[iAlpha] = lookupCd.find(alphaNew[iAlpha])

    # Calculate the initial step's variables (using the current gamma and
    # atmospheric properties)
    vVerNew = vVerCur + 0.5 * ( K1 * (
        (ClCosCur + ClNew * np.cos(gammaCur)) +
        (CdSinCur + CdNew * np.sin(gammaCur)) )
        - 2 * gCur
    ) * dt

    vHorNew = vHorCur + 0.5 * K1 * (
        (ClSinCur + ClNew * np.sin(gammaCur)) -
        (CdCosCur + CdNew * np.cos(gammaCur))
    ) * dt

    # Start iterating until a stable solution is found
    for iIteration in range(0, 100):
        # Calculate the new gamma and atmospheric properties
        hNew = hCur + 0.5 * (vVerCur + vVerNew) * dt
        vWindZonalNew = atmosphere.velocityZonal(hNew, latitudeCur, longitudeCur)[1]
        gammaNew = np.arctan2(vVerNew, vWindZonalNew + vHorNew)
        rhoNew = atmosphere.density(hNew, latitudeCur, longitudeCur)[1]
        vInfNewSquared = np.power(vHorNew + vWindZonalNew, 2.0) + np.power(vVerNew, 2.0)
        K2 = gNew * 0.5 * rhoNew * vInfNewSquared * S / W

        # Calculate the new velocities with the updated angle and atmospheric
        # estimations
        vVerNewGuess = vVerCur + 0.5 * ((
            K1 * (ClCosCur + CdCosCur) +
            K2 * (ClNew * np.cos(gammaNew) + CdNew * np.sin(gammaNew))
        ) - (gCur + gNew)) * dt

        vHorNewGuess = vHorCur + 0.5 * (
            K1 * (ClSinCur - CdCosCur) +
            K2 * (ClNew * np.sin(gammaNew) - CdNew * np.cos(gammaNew))
        ) * dt

        # Check if the new estimations are all below the desired threshold
        vFactorVer = (vVerNewGuess - vVerNew) / vVerNew
        converged = True

        for iFactor in range(0, len(vFactorVer)):
            if vFactorVer[iFactor] > tol:
                # Vertical speed has not converged enough
                #print('vFactorVer', iFactor, '=', vFactorVer[iFactor], 'too large ( iteration =', iIteration, ')')
                converged = False
                break

        if not converged:
            vHorNew = vHorNewGuess
            vVerNew = vVerNewGuess
            continue
        
        vFactorHor = (vHorNewGuess - vHorNew) / vHorNew

        for iFactor in range(0, len(vFactorHor)):
            if vFactorHor[iFactor] > tol:
                # Horizontal speed has not converged enough
                #print('vFactorHor', iFactor, '=', vFactorHor[iFactor], 'too large (iteration =', iIteration, ')')
                converged = False
                break
            
        if not converged:
            vHorNew = vHorNewGuess
            vVerNew = vVerNewGuess
            break

        # If this position is reached then all values have converged to the
        # desired degree
        #print('Solution converged after', iIteration + 1, 'iterations')
        return vHorNewGuess, vVerNewGuess, gammaNew

    raise RuntimeError("TrackDive.Step did not converge")

# Analysis_GenerateLiftOverDrag generates two lookup tables to lookup the
# aerodynamic performance parameter L/D. The first item in the returned is the
# forward lookup table (use an angle of attack to lookup a L/D value). The
# second item is the backward lookup table (use a L/D value to lookup an angle
# of attack). The second item will return arrays of values (and can return
# two angels of attack). The third and the fourth item contain lookup tables
# for the lift coefficient (forward and backward)
def TestGenerateLiftOverDrag():
    dataAlpha = np.linspace(-5.0, 5.0, 21)
    dataClCd = [
        -24.0, -23.0, -23.3, -22.1, -20.1, -17.3, -13.5, -8.8, -3.3, 2.8, 8.8,
        14.3, 18.9, 22.3, 24.6, 26.1, 26.6, 26.6, 26.1, 25.5, 24.7
    ]
    dataCl = [
        -0.333, -0.294, -0.255, -0.216, -0.177, -0.139, -0.100, -0.0605,
        -0.0214, 0.0177, 0.0568, 0.0959, 0.134, 0.174, 0.213, 0.252,
        0.291, 0.329, 0.368, 0.407, 0.445
    ]
    dataCd = [
        0.0138, 0.0123, 0.0110, 0.00980, 0.00881, 0.00800, 0.00738,
        0.00689, 0.00657, 0.00638, 0.00642, 0.00668, 0.00713, 0.00778,
        0.00865, 0.00966, 0.0109, 0.0124, 0.0141, 0.0160, 0.0180
    ]

    lutPerfForward = TrackLookup.Lookup1D(dataAlpha, dataClCd)
    lutPerfBackward = TrackLookup.LookupSegmented1D(dataClCd, dataAlpha)

    lutClForward = TrackLookup.Lookup1D(dataAlpha, dataCl)
    lutClBackward = TrackLookup.LookupSegmented1D(dataCl, dataAlpha)

    lutCdForward = TrackLookup.Lookup1D(dataAlpha, dataCd)
    lutCdBackward = TrackLookup.LookupSegmented1D(dataCd, dataAlpha)

    return lutPerfForward, lutPerfBackward, \
        lutClForward, lutClBackward, \
        lutCdForward, lutCdBackward

# TestStep is the testing diving function that will plot the results for a
# given set of variables.
def TestStep(h, alpha, gamma, alphaNew, Vhor, Vver, dt, W, S, luts, atm):
    # Calculate often used variables
    rho = atm.density(h, 0, 0)
    vZonal = atm.velocityZonal(h, 0, 0)
    g1 = 8.8
    g2 = 8.8
    Vinf = np.sqrt((Vhor + vZonal[1])**2.0 + Vver**2.0)
    K1 = g1 * 0.5 * rho[1] * Vinf**2.0 * S / W
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

    # Calculate the initial step
    VverNew = Vver + 0.5 * ( K1 * (
        (ClCos1 + Cl2 * np.cos(gamma)) +
        (CdSin1 + Cd2 * np.sin(gamma)) )
        - 2 * g1
    ) * dt

    VhorNew = Vhor + 0.5 * K1 * (
        (ClSin1 + Cl2 * np.sin(gamma)) -
        (CdCos1 + Cd2 * np.cos(gamma))
    ) * dt

    dhNew = 0.5 * (Vver + VverNew) * dt

    vZonalNew = atm.velocityZonal(h + dhNew, 0, 0)

    gammaNew = np.arctan2(VverNew, vZonalNew[1] + VhorNew)

    # Start looping to find a stable solution
    K2 = np.zeros(alphaNew.shape)

    VverAll = []
    VhorAll = []
    VverAll.append(VverNew)
    VhorAll.append(VhorNew)

    for i in range(0, 10):
        print('ITERATION ---------', i, ', dt =', dt)
        hNew = h + dhNew
        rhoNew = atm.density(hNew, 0, 0)
        VinfNew = np.sqrt(np.power(VhorNew + vZonalNew[1], 2.0) + np.power(VverNew, 2.0))
        K2 = g2 * 0.5 * rhoNew[1] * np.power(VinfNew, 2.0) * S / W

        print('gammaNew =', gammaNew)
        print('term1 =', K1 * (ClCos1 + CdSin1))
        print('term2 =', K2 * (Cl2 * np.cos(gammaNew) + Cd2 * np.sin(gammaNew)))
        print('term3 =', (g1 + g2))

        VverNew = Vver + 0.5 * ((
            K1 * (ClCos1 + CdSin1) +
            K2 * (Cl2 * np.cos(gammaNew) + Cd2 * np.sin(gammaNew))
        ) - (g1 + g2)) * dt

        VhorNew = Vhor + 0.5 * (
            K1 * (ClSin1 - CdCos1) +
            K2 * (Cl2 * np.sin(gammaNew) - Cd2 * np.cos(gammaNew))
        ) * dt

        dhNew = 0.5 * (Vver + VverNew) * dt
        print('Vver =', Vver)
        print('VverNew =', VverNew)
        print('h =', h)
        print('dhNew =', dhNew)
        vZonalNew = atm.velocityZonal(h + dhNew, 0, 0)
        gammaNew = np.arctan2(VverNew, vZonalNew[1] + VhorNew)

        VverAll.append(VverNew)
        VhorAll.append(VhorNew)

    VverAll = np.asarray(VverAll)
    VhorAll = np.asarray(VhorAll)

    fig = plt.figure()
    axVer = fig.add_subplot(121)
    axHor = fig.add_subplot(122)
    linesVer = []
    linesHor = []
    name = []
    cmap = plt.get_cmap('jet')

    for i in range(0, VverAll.shape[0]):
        color = cmap(i / (VverAll.shape[0] - 1))
        lVer, = axVer.plot(alphaNew, VverAll[i] - VverAll[-1], color=color)
        lHor, = axHor.plot(alphaNew, VhorAll[i] - VhorAll[-1], color=color)

        linesVer.append(lVer)
        linesHor.append(lHor)
        name.append('i = ' + str(i + 1))

    axVer.legend(linesVer, name)
    axHor.legend(linesHor, name)
    axVer.set_title('vertical')
    axHor.set_title('horizontal')

#atm = Atmosphere.Atmosphere()
#luts = TestGenerateLiftOverDrag()
#TestStep(55000, -2.0, -5.0 / 180.0 * np.pi, np.linspace(-4.0, 4.0, 10), 50, -5, 2, 1, 0.001, luts, atm)
#Step(55000, -2.0, -5.0 / 180.0 * np.pi, 30, -5, 0, 0, 1, 0.001, np.linspace(-5.0, 5.0, 10),
#    1.0, luts[2], luts[4], atm, 1e-15)
