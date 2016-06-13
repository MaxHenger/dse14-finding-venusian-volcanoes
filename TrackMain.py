#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 11:08:41 2016

@author: MaxHenger
"""

import TrackBounds
import TrackSettings
import TrackCommon
import TrackStitching
import TrackClimbOptimize
import TrackLookup
import Atmosphere

import numpy as np

def Main():
    # Determine track operational points
    axisHeight = np.linspace(30e3, 80e3, 50)
    axisDeltaV = np.linspace(-50, 50, 50)
    axisSeverity = np.linspace(-2.5, 2.5, 50)
    
    settings = TrackSettings.Settings()
    atmosphere = Atmosphere.Atmosphere()
    PReqMin = 10
    PReq = 32e3
    dt = 0.25

    operationalRegions = TrackBounds.DetermineTracks(axisHeight, axisDeltaV,
        settings.latitude, settings.longitude, PReq, axisSeverity)
    
    print(' > Found', len(operationalRegions), 'operational regions')
    
    middle = int(len(operationalRegions) / 2)
    operationalRegions = [operationalRegions[middle]]
    
    # For each of the operational regions simulate the track
    severity = []
    upperHeight = []
    upperDeltaV = []
    lowerHeight = []
    lowerDeltaV = []
    result = []
    simulationFilename = []

    for iOperationalRegion in range(0, len(operationalRegions)):
        # For each of the lower and upper operational points simulate the 
        # response.
        curOperationalRegion = operationalRegions[iOperationalRegion]
        curSeverity = curOperationalRegion[0]
        curUpperDeltaV = curOperationalRegion[1]
        curUpperHeight = curOperationalRegion[2]
        curLowerDeltaV = curOperationalRegion[3]
        curLowerHeight = curOperationalRegion[4]

        # Generate ascent guide file
        print(TrackCommon.StringHeader("Determining ascent guides", 60))
        print(TrackCommon.StringPad(" > Severity:", curSeverity, 2, 6))
        ascentGuide = TrackClimbOptimize.GenerateAscentMaps(axisHeight, 
            axisDeltaV, settings.W, settings.S, settings.inclination,
            settings.lookupCl, settings.lookupCd, atmosphere, settings.qInfMin,
            settings.qInfMax, settings.alphaMin, settings.alphaMax, PReqMin, PReq,
            severity=curSeverity)
        
        lowerGuide = TrackLookup.Lookup1D(ascentGuide['pathHeight'], 
            ascentGuide['pathMinDeltaV'])

        settings.lowerBound = lowerGuide
        settings.upperBound = None
        
        # Start iterating through all operational points
        for iUpper in range(0, 2):
            for iLower in range(0, 1):
        #for iUpper in range(0, len(curUpperDeltaV)):
        #    for iLower in range(0, len(curLowerDeltaV)):
                # Print information
                print(TrackCommon.StringHeader("Simulating response", 60))
                print(TrackCommon.StringPad("severity:    ", curSeverity, 3, 6))
                print(TrackCommon.StringPad("upper height:", curUpperHeight[iUpper] / 1e3, 2, 6), "km")
                print(TrackCommon.StringPad("upper deltaV:", curUpperDeltaV[iUpper], 2, 6), "m/s")
                print(TrackCommon.StringPad("lower height:", curLowerHeight[iLower] / 1e3, 2, 6), "km")
                print(TrackCommon.StringPad("lower deltaV:", curLowerDeltaV[iLower], 2, 6), "m/s")
                
                filename = ''
                success = False
                
                try:
                    # Attempt to simulate track. It will throw an exception if
                    # somehow the solution cannot be simulated or will not
                    # converge
                    filename = TrackStitching.StitchTracks(
                        curUpperHeight[iUpper], curUpperDeltaV[iUpper], 
                        curLowerHeight[iLower], curLowerDeltaV[iLower], 10, 
                        curLowerDeltaV[iLower], curUpperDeltaV[iUpper], PReqMin, 
                        PReq, dt, settings, curSeverity, plotResult=False, 
                        saveResult=True)
                    
                    success = True
                except Exception as ex:
                    print(TrackCommon.StringHeader("ERROR: " + str(ex), 60, '!', '!', '!'))
                    success = False
                    
                severity.append(curSeverity)
                upperHeight.append(curUpperHeight[iUpper])
                upperDeltaV.append(curUpperDeltaV[iUpper])
                lowerHeight.append(curLowerHeight[iLower])
                lowerDeltaV.append(curLowerDeltaV[iLower])
                
                if success:
                    result.append(1)
                else:
                    result.append(0)
                    
                simulationFilename.append(filename)
    
    # Save all final results to a csv file
    np.savetxt('final.csv', list(map(list, zip(*[severity, upperHeight, 
        upperDeltaV, lowerHeight, lowerDeltaV, result, simulationFilename]))),
        ['%.10f', '%.10f', '%.10f', '%.10f', '%.10f', '%d', '%s'],
        delimiter=';', header='severity; upper height [m]; upper deltaV [m/s]; ' + 
            'lower height [m]; lower deltaV [m/s]; result [1/0]; filename')
                        
Main()