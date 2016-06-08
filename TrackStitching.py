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
import TrackSettings
import TrackStorage

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as scp_int

def GroundRotation(time, height, vVer, settings):
    y = vVer / (settings.RVenus + height) - settings.omegaVenus
    return scp_int.trapz(y, time)

def SetAxisColors(ax, color):
    for tick in ax.get_yticklabels():
        tick.set_color(color)

def StitchTracks(preDiveHeight, preDiveVHor, postDiveHeight, postDiveVHor,
                 postDiveLoiter, preAscentVHor, postAscentVHor,
                 PReqMin, PReqMax, dt, settings, severity, saveResult=True):
    # For less verbose typing
    W = settings.W
    S = settings.S
    latitude = settings.latitude
    longitude = settings.longitude
    inclination = settings.inclination

    # Load all required data
    atm = Atmosphere.Atmosphere()
    lookupCl, lookupCd = TrackCommon.LoadAerodynamicData(
        './data/aerodynamicPerformance/Cl.csv', './data/aerodynamicPerformance/Cd.csv')
    lookupLowerAscent, lookupUpperAscent = TrackCommon.LoadAscentGuides(
        'optclimb_-60.0to20.0_0.0.dat', 2.0)

    # Start by performing a dive
    timeDive, heightDive, vHorDive, vVerDive, vInfDive, alphaDive, gammaDive = \
        TrackDiveOptimize.OptimizeDive(preDiveHeight, postDiveHeight,
        preDiveVHor, 0, longitude, latitude, W, S, postDiveVHor, 0, dt,
        lookupCl, lookupCd, severity, plotResults=False, storeResults=False)
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
        plotResults=False, storeResults=False)
    heightAcc1 = np.repeat(heightDive[-1], len(timeAcc1))
    vVerAcc1 = np.zeros([len(timeAcc1)])
    vInfAcc1 = vZonalAcc1 + vHorAcc1
    gammaAcc1 = np.zeros([len(timeAcc1)])

    timeEndAcc1 = timeAcc1[-1] + timeEndDive + dt

    # Performing loiter. Appending these arrays (inefficiently) because one day
    # I might actually implement some proper loitering
    postDiveLoiterNum = int(postDiveLoiter / dt)
    if postDiveLoiterNum == 0:
        # waste at least one dt to make writing code easier
        postDiveLoiterNum = 1

    timeLoiter1 = np.linspace(0, (postDiveLoiterNum - 1) * dt, postDiveLoiterNum)
    heightLoiter1 = np.repeat(heightAcc1[-1], postDiveLoiterNum)
    vHorLoiter1 = np.repeat(vHorAcc1[-1], postDiveLoiterNum)
    vVerLoiter1 = np.zeros([postDiveLoiterNum])
    vInfLoiter1 = np.repeat(vInfAcc1[-1], postDiveLoiterNum)
    alphaLoiter1 = np.repeat(alphaAcc1[-1], postDiveLoiterNum)
    gammaLoiter1 = np.repeat(gammaAcc1[-1], postDiveLoiterNum)
    powerLoiter1 = np.repeat(powerAcc1[-1], postDiveLoiterNum)

    timeEndLoiter1 = timeLoiter1[-1] + timeEndAcc1 + dt

    # Perform the second speed change before initiating the climb
    vZonalAcc2 = TrackCommon.AdjustSeverity(atm.density(heightLoiter1[-1],
        latitude, longitude), severity)

    timeAcc2, vHorAcc2, alphaAcc2, powerAcc2 = \
        TrackAcceleratingOptimize.OptimizeAccelerating(heightLoiter1[-1],
        vHorLoiter1[-1], powerLoiter1[-1], alphaLoiter1[-1], longitude,
        latitude, W, S, preAscentVHor, inclination, dt, PReqMin, PReqMax,
        lookupCl, lookupCd, severity, plotResults=False, storeResults=False)
    heightAcc2 = np.repeat(heightLoiter1[-1], len(timeAcc2))
    vVerAcc2 = np.zeros([len(timeAcc2)])
    vInfAcc2 = vZonalAcc2 + vHorAcc2
    gammaAcc2 = np.zeros([len(timeAcc2)])

    timeEndAcc2 = timeAcc2[-1] + timeEndLoiter1 + dt

    # Start the ascent
    heightClimbQuit = heightAcc1[-1] - 0.1 * (preDiveHeight - heightAcc1[-1])

    timeClimb, heightClimb, vHorClimb, vVerClimb, vInfClimb, alphaClimb, gammaClimb = \
        TrackClimbOptimize.OptimizeClimb(heightAcc2[-1], preDiveHeight,
        heightClimbQuit, vHorAcc2[-1], vVerAcc2[-1], longitude, latitude, W, S,
        postAscentVHor, 0, PReqMax, inclination, dt, lookupCl, lookupCd,
        severity, lookupBoundLowerVInf=lookupLowerAscent,
        lookupBoundUpperVInf=lookupUpperAscent, plotResults=False,
        storeResults=False)
    powerClimb = np.repeat(PReqMax, len(timeClimb))

    timeEndClimb = timeClimb[-1] + timeEndAcc2 + dt

    # Accelerate to post ascent horizontal velocity
    vZonalAcc3 = TrackCommon.AdjustSeverity(atm.density(heightClimb[-1],
        latitude, longitude), severity)

    timeAcc3, vHorAcc3, alphaAcc3, powerAcc3 = \
        TrackAcceleratingOptimize.OptimizeAccelerating(heightClimb[-1],
        vHorClimb[-1], powerClimb[-1], alphaClimb[-1], longitude, latitude,
        W, S, postAscentVHor, inclination, dt, PReqMin, PReqMax, lookupCl, lookupCd,
        severity, plotResults=False, storeResults=False)
    heightAcc3 = np.repeat(heightClimb[-1], len(timeAcc3))
    vVerAcc3 = np.zeros([len(timeAcc3)])
    vInfAcc3 = vZonalAcc3 + vHorAcc3
    gammaAcc3 = np.zeros([len(timeAcc3)])

    timeEndAcc3 = timeAcc3[-1] + timeEndClimb + dt

    # Already perform the post-loiter acceleration before diving. This is used
    # to estimate the time needed at loitering to end up at the same subsolar
    # point in the end
    vZonalAcc4 = TrackCommon.AdjustSeverity(atm.density(heightAcc3[-1],
        latitude, longitude), severity)

    timeAcc4, vHorAcc4, alphaAcc4, powerAcc4 = \
        TrackAcceleratingOptimize.OptimizeAccelerating(heightAcc3[-1],
        vHorAcc3[-1], powerAcc3[-1], alphaAcc3[-1], longitude, latitude,
        W, S, preDiveVHor, inclination, dt, PReqMin, PReqMax, lookupCl, lookupCd,
        severity, plotResults=False, storeResults=False)
    heightAcc4 = np.repeat(heightAcc3[-1], len(timeAcc4))
    vVerAcc4 = np.zeros([len(timeAcc4)])
    vInfAcc4 = vZonalAcc4 + vHorAcc4
    gammaAcc4 = np.zeros([len(timeAcc4)])

    # Determine the upper-height loiter time required to reach the same subsolar
    # point again.
    rotationCovered = 0.0

    for t, h, v in [(timeDive, heightDive, vHorDive),
                    (timeAcc1, heightAcc1, vHorAcc1),
                    (timeLoiter1, heightLoiter1, vHorLoiter1),
                    (timeAcc2, heightAcc2, vHorAcc2),
                    (timeClimb, heightClimb, vHorClimb),
                    (timeAcc3, heightAcc3, vHorAcc3),
                    (timeAcc4, heightAcc4, vHorAcc4)]:
        rotationCovered += GroundRotation(np.asarray(t), np.asarray(h),
                                          np.asarray(v), settings)

    timeLoiterUp = ((settings.RVenus + heightAcc3[-1]) / (settings.omegaVenus *
        (settings.RVenus + heightAcc3[-1]) - vHorAcc3[-1]) * rotationCovered)

    if timeLoiterUp < 0:
        raise ValueError("Upper loiter time is smaller than 0:", timeLoiterUp)

    postClimbLoiterNum = int(timeLoiterUp / dt)

    print(' > spending', round(timeLoiterUp, 1), 's at top (', postClimbLoiterNum, 'iterations )')

    if postClimbLoiterNum == 0:
        postClimbLoiterNum = 1

    timeLoiter2 = np.linspace(0, (postClimbLoiterNum - 1) * dt, postClimbLoiterNum)
    heightLoiter2 = np.repeat(heightAcc3[-1], postClimbLoiterNum)
    vHorLoiter2 = np.repeat(vHorAcc3[-1], postClimbLoiterNum)
    vVerLoiter2 = np.zeros([postClimbLoiterNum])
    vInfLoiter2 = np.repeat(vInfAcc3[-1], postClimbLoiterNum)
    alphaLoiter2 = np.repeat(alphaAcc3[-1], postClimbLoiterNum)
    gammaLoiter2 = np.repeat(gammaAcc3[-1], postClimbLoiterNum)
    powerLoiter2 = np.repeat(powerAcc3[-1], postClimbLoiterNum)

    timeEndLoiter2 = timeLoiter2[-1] + timeEndAcc3 + dt
    timeEndAcc4 = timeAcc4[-1] + timeEndLoiter2 + dt

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

    timeTotal.extend(np.asarray(timeAcc3) + timeEndClimb)
    heightTotal.extend(heightAcc3)
    vHorTotal.extend(vHorAcc3)
    vVerTotal.extend(vVerAcc3)
    vInfTotal.extend(vInfAcc3)
    alphaTotal.extend(alphaAcc3)
    gammaTotal.extend(gammaAcc3)

    timeTotal.extend(np.asarray(timeLoiter2) + timeEndAcc3)
    heightTotal.extend(heightLoiter2)
    vHorTotal.extend(vHorLoiter2)
    vVerTotal.extend(vVerLoiter2)
    vInfTotal.extend(vInfLoiter2)
    alphaTotal.extend(alphaLoiter2)
    gammaTotal.extend(gammaLoiter2)

    timeTotal.extend(np.asarray(timeAcc4) + timeEndLoiter2)
    heightTotal.extend(heightAcc4)
    vHorTotal.extend(vHorAcc4)
    vVerTotal.extend(vVerAcc4)
    vInfTotal.extend(vInfAcc4)
    alphaTotal.extend(alphaAcc4)
    gammaTotal.extend(gammaAcc4)

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

    fig.suptitle('Flight parameters')

    # Store all data
    if saveResult:
        file = TrackStorage.DataStorage()
        file.addVariable('timeTotal', timeTotal)
        file.addVariable('heightTotal', heightTotal)
        file.addVariable('vHorTotal', vHorTotal)
        file.addVariable('vVerTotal', vVerTotal)
        file.addVariable('vInfTotal', vInfTotal)
        file.addVariable('alphaTotal', alphaTotal)
        file.addVariable('gammaTotal', gammaTotal)
        file.save('stitched_' + str(round(preDiveHeight, 1)) + "at" +
                  str(round(preDiveVHor, 1)) + 'to' + str(round(postDiveHeight, 1)) +
                  'at' + str(round(postDiveHeight, 1)) + '.dat')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    groundCovered = TrackCommon.CumulativeSimps(vHorTotal * settings.RVenus /
        (settings.RVenus + heightTotal) - settings.omegaVenus * settings.RVenus,
        timeTotal)

    ax.plot(groundCovered, heightTotal)

StitchTracks(62e3, 3.5, 38e3, -25.0, 10, -5.0, 7.8,
             0e3, 32e3, 0.20, TrackSettings.Settings(), 0.0)
