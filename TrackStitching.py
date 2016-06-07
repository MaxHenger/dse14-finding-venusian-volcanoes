# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 00:58:02 2016

@author: MaxHenger
"""

import Atmosphere
import TrackCommon
import TrackDiveOptimize
import TrackAcceleratingOptimize
import TrackClimbOptimize

import matplotlib.pyplot as plt
import numpy as np

def SetAxisColors(ax, color):
    for tick in ax.get_yticklabels():
        tick.set_color(color)

def StitchTracks(preDiveHeight, preDiveVHor, postDiveHeight, postDiveVHor,
                 postDiveLoiter, preAscentVHor, postAscentVHor, postAscentLoiter,
                 W, S, latitude, longitude, inclination, PReqMin, PReqMax, dt,
                 severity):
    # Load all required data
    atm = Atmosphere.Atmosphere()
    lookupCl, lookupCd = TrackCommon.LoadAerodynamicData(
        './data/aerodynamicPerformance/Cl.csv', './data/aerodynamicPerformance/Cd.csv')
    lookupLowerAscent, lookupUpperAscent = TrackCommon.LoadAscentGuides(
        'optclimb_-60.0to20.0_0.0.dat')

    # Start by performing a dive
    timeDive, heightDive, vHorDive, vVerDive, vInfDive, alphaDive, gammaDive = \
        TrackDiveOptimize.OptimizeDive(preDiveHeight, postDiveHeight,
        preDiveVHor, 0, longitude, latitude, W, S, postDiveVHor, 0, dt,
        lookupCl, lookupCd, severity, plotResults=True, storeResults=True)
    powerDive = np.zeros([len(timeDive)])

    timeEndDive = timeDive[-1] + dt

    # Change the speed to the desired one, use the final values from the dive
    # as input
    vZonalAcc1 = TrackCommon.AdjustSeverity(atm.density(heightDive[-1],
        latitude, longitude), severity)

    timeAcc1, vHorAcc1, alphaAcc1, powerAcc1 = \
        TrackAcceleratingOptimize.OptimizeAccelerating(heightDive[-1],
        vHorDive[-1], 0, alphaDive[-1], longitude, latitude, W, S, postDiveVHor,
        inclination, dt, PReqMin, PReqMax, lookupCl, lookupCd, severity,
        plotResults=True, storeResults=True)
    heightAcc1 = np.repeat(heightDive[-1], len(timeAcc1))
    vVerAcc1 = np.zeros([len(timeAcc1)])
    vInfAcc1 = vZonalAcc1 + vHorAcc1
    gammaAcc1 = np.zeros([len(timeAcc1)])

    timeEndAcc1 = timeAcc1[-1] + timeEndDive + dt

    # Performing loiter. Appending these arrays (inefficiently) because one day
    # I might actually implement some proper loitering
    postAscentLoiterNum = int(postAscentLoiter / dt)
    if postAscentLoiterNum == 0:
        # waste at least one deltaV to make writing code easier
        postAscentLoiterNum = 1
    timeLoiter1 = np.linspace(0, (postAscentLoiterNum - 1) * dt, postAscentLoiterNum)
    heightLoiter1 = np.repeat(heightAcc1[-1], postAscentLoiterNum)
    vHorLoiter1 = np.repeat(vHorAcc1[-1], postAscentLoiterNum)
    vVerLoiter1 = np.zeros([postAscentLoiterNum])
    vInfLoiter1 = np.repeat(vInfAcc1[-1], postAscentLoiterNum)
    alphaLoiter1 = np.repeat(alphaAcc1[-1], postAscentLoiterNum)
    gammaLoiter1 = np.repeat(gammaAcc1[-1], postAscentLoiterNum)
    powerLoiter1 = np.repeat(powerAcc1[-1], postAscentLoiterNum)

    timeEndLoiter1 = timeLoiter1[-1] + timeEndAcc1 + dt

    # Perform the second speed change before initiating the climb
    vZonalAcc2 = TrackCommon.AdjustSeverity(atm.density(heightLoiter1[-1],
        latitude, longitude), severity)

    timeAcc2, vHorAcc2, alphaAcc2, powerAcc2 = \
        TrackAcceleratingOptimize.OptimizeAccelerating(heightLoiter1[-1],
        vHorLoiter1[-1], powerLoiter1[-1], alphaLoiter1[-1], longitude,
        latitude, W, S, preAscentVHor, inclination, dt, PReqMin, PReqMax, 
        lookupCl, lookupCd, severity, plotResults=True, storeResults=True)
    heightAcc2 = np.repeat(heightLoiter1[-1], len(timeAcc2))
    vVerAcc2 = np.zeros([len(timeAcc2)])
    vInfAcc2 = vZonalAcc2 + vHorAcc2
    gammaAcc2 = np.zeros([len(timeAcc2)])

    timeEndAcc2 = timeAcc2[-1] + timeEndLoiter1 + dt

    # Start the ascent
    heightClimbQuit = heightAcc1[-1] - 0.1 * (preDiveHeight - heightAcc1[-1])
    timeClimb, heightClimb, vHorClimb, vVerClimb, vInfClimb, alphaClimb, gammaClimb = \
        TrackClimbOptimize.OptimizeClimb(heightAcc1[-1], preDiveHeight,
        heightClimbQuit, vHorAcc2[-1], vVerAcc2[-1], longitude, latitude, W, S,
        postAscentVHor, 0, PReqMax, inclination, dt, lookupCl, lookupCd,
        severity, lookupBoundLowerVInf=lookupLowerAscent,
        lookupBoundUpperVInf=lookupUpperAscent, plotResults=True,
        storeResults=True)

    timeEndClimb = timeClimb[-1] + timeEndAcc2 + dt

    # Loiter for as long as planned by simply repeating the post-acceleration
    # values for as long as required

    # Compound all data
    # - create empty arrays
    timeTotal = []
    heightTotal = []
    vHorTotal = []
    vVerTotal = []
    vInfTotal = []
    alphaTotal = []
    gammaTotal = []

    # - glue all data together in a plottable manner
    timeTotal.extend(timeDive)
    heightTotal.extend(heightDive)
    vHorTotal.extend(vHorDive)
    vVerTotal.extend(vVerDive)
    vInfTotal.extend(vInfDive)
    alphaTotal.extend(alphaDive)
    gammaTotal.extend(gammaDive)

    timeTotal.extend(np.asarray(timeAcc1) + timeEndDive)
    heightTotal.extend(heightAcc1)
    vHorTotal.extend(vHorAcc1)
    vVerTotal.extend(vVerAcc1)
    vInfTotal.extend(vInfAcc1)
    alphaTotal.extend(alphaAcc1)
    gammaTotal.extend(gammaAcc1)

    timeTotal.extend(np.asarray(timeLoiter1) + timeEndAcc1)
    heightTotal.extend(heightLoiter1)
    vHorTotal.extend(vHorLoiter1)
    vVerTotal.extend(vVerLoiter1)
    vInfTotal.extend(vInfLoiter1)
    alphaTotal.extend(alphaLoiter1)
    gammaTotal.extend(gammaLoiter1)

    timeTotal.extend(np.asarray(timeAcc2) + timeEndLoiter1)
    heightTotal.extend(heightAcc2)
    vHorTotal.extend(vHorAcc2)
    vVerTotal.extend(vVerAcc2)
    vInfTotal.extend(vInfAcc2)
    alphaTotal.extend(alphaAcc2)
    gammaTotal.extend(gammaAcc2)

    timeTotal.extend(np.asarray(timeClimb) + timeEndAcc2)
    heightTotal.extend(heightClimb)
    vHorTotal.extend(vHorClimb)
    vVerTotal.extend(vVerClimb)
    vInfTotal.extend(vInfClimb)
    alphaTotal.extend(alphaClimb)
    gammaTotal.extend(gammaClimb)

    # - convert to numpy arrays for data manipulation
    timeTotal = np.asarray(timeTotal)
    heightTotal = np.asarray(heightTotal)
    vHorTotal = np.asarray(vHorTotal)
    vVerTotal = np.asarray(vVerTotal)
    vInfTotal = np.asarray(vInfTotal)
    alphaTotal = np.asarray(alphaTotal)
    gammaTotal = np.asarray(gammaTotal)

    # Start plotting
    fig = plt.figure()
    axHeight = fig.add_subplot(311)
    axVelCom = fig.add_subplot(312)
    axVelInf = axVelCom.twinx()
    axAlpha = fig.add_subplot(313)
    axGamma = axAlpha.twinx()

    axHeight.plot(timeTotal, heightTotal / 1e3, 'g')
    axHeight.set_xlabel('time [s]')
    axHeight.set_ylabel('height [km]')
    axHeight.grid(True)

    axVelCom.plot(timeTotal, vHorTotal, 'r', label='vHor')
    axVelCom.plot(timeTotal, vVerTotal, 'r--', label='vVer')
    axVelInf.plot(timeTotal, vInfTotal, 'g', label='vInf')

    SetAxisColors(axVelCom, 'r')
    SetAxisColors(axVelInf, 'g')

    axVelCom.set_xlabel('time [s]')
    axVelCom.set_ylabel('velocity [m/s]')
    axVelInf.set_ylabel('velocity [m/s]')
    axVelCom.legend()
    axVelCom.grid(True)

    axAlpha.plot(timeTotal, alphaTotal, 'r', label='alpha')
    axGamma.plot(timeTotal, gammaTotal, 'g', label='gamma')

    SetAxisColors(axAlpha, 'r')
    SetAxisColors(axGamma, 'g')

    axAlpha.set_xlabel('time [s]')
    axAlpha.set_ylabel('alpha [deg]')
    axGamma.set_ylabel('gamma [deg]')

StitchTracks(62e3, 3.5, 32e3, -25.0, 0,
             -25.0, 7.8, 0,
             700*8.8, 35, 0, 0, 0, 0e3, 32e3, 0.35,
             0.0)
