#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 02:45:19 2016

@author: MaxHenger
"""

import Atmosphere
import TrackCommon

import numpy as np
import matplotlib.pyplot as plt

def PlotQInfty():
    atm = Atmosphere.Atmosphere()
    deltaV = np.linspace(-60, 60, 241)
    height = np.linspace(30e3, 80e3, 501)
    qInftyMin = np.zeros([len(height), len(deltaV)])
    qInfty = [qInftyMin, np.zeros(qInftyMin.shape), np.zeros(qInftyMin.shape)]
    severities = [-1.5, 0, 1.5]

    for iSeverity in range(0, len(severities)):
        density = TrackCommon.AdjustSeverity(atm.density(
            height, 0, 0), severities[iSeverity])
        vZonal = TrackCommon.AdjustSeverity(atm.velocityZonal(
            height, 0, 0), severities[iSeverity])
        
        for iDeltaV in range(0, len(deltaV)):
            vInf = vZonal + deltaV[iDeltaV]
            qInf = 0.5 * density * np.power(vInf, 2.0)
            
            qInfty[iSeverity][:, iDeltaV] = np.log10(qInf)
        
        boundsData, boundsLegend = TrackCommon.ImageAxes(0.07, 0.93, 0.04, 1.0)
        fig = plt.figure()
        axData = fig.add_axes(boundsData)
        axLegend = fig.add_axes(boundsLegend)
        TrackCommon.PlotImage(fig, axData, axLegend, deltaV, r'$V_I \; [m/s]$',
            height / 1e3, r'$h \; [km]$', qInfty[iSeverity], 
            r'$\mathrm{log}_{10} ( q_\infty ) \; [Pa]$',
            contours=[1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0], forceNormMin=0.1)
        axData.grid(True)
        
PlotQInfty()