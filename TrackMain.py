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
import TrackStorage
import TrackPower
import AircraftBatteries
import Atmosphere
import TimeEstimator

import numpy as np
import os.path as os_p
import scipy.integrate as scp_int
import matplotlib.pyplot as plt

def Main():
    filename = 'final.csv'
    RESULT_PENDING = 0
    RESULT_DONE = 1
    RESULT_ERROR = 2

    # Determine track operational points
    axisHeight = np.linspace(30e3, 80e3, 75)
    axisDeltaV = np.linspace(-100, 50, 100)
    axisSeverity = np.linspace(-2.5, 2.5, 75)

    settings = TrackSettings.Settings()
    atmosphere = Atmosphere.Atmosphere()
    powerEstimator = TrackPower.TrackPower(settings, atmosphere)

    PReqMin = 10
    PReq = 32e3
    dt = 0.25
    climbLowerBoundOffset = 0.0

    # Check if an operational regions file was created before
    severity = []
    upperHeight = []
    upperDeltaV = []
    lowerHeight = []
    lowerDeltaV = []
    result = []
    simulationFilename = []
    areaTop = []
    batteryCapacity = []
    batteryWeight = []

    toProcess = 0
    
    fileHeader = 'severity; upper height [m]; upper deltaV [m/s]; lower height [m]; ' + \
                 'lower deltaV [m/s]; success; filename; areaTop [m2]; ' + \
                 'capacity [J]; battery mass [kg]\n'

    if os_p.isfile(filename):
        # Read csv file
        fh = open(filename, 'r')
        line = fh.readline() # first line is the header
        line = fh.readline()

        while len(line) != 0:
            # Process the current line
            seperated = line.split(';')

            for iSeperated in range(0, len(seperated)):
                seperated[iSeperated] = seperated[iSeperated].strip()

            severity.append(float(seperated[0]))
            upperHeight.append(float(seperated[1]))
            upperDeltaV.append(float(seperated[2]))
            lowerHeight.append(float(seperated[3]))
            lowerDeltaV.append(float(seperated[4]))
            result.append(int(seperated[5]))
            simulationFilename.append(seperated[6])
            areaTop.append(seperated[7])
            batteryCapacity.append(seperated[8])
            batteryWeight.append(seperated[9])

            if result[-1] == RESULT_PENDING:
                toProcess += 1

            line = fh.readline()

        fh.close()
    else:
        # No file existed yet
        operationalRegions = TrackBounds.DetermineTracks(axisHeight, axisDeltaV,
            settings.latitude, settings.longitude, PReq, axisSeverity)

        print(' > Found', len(operationalRegions), 'operational regions')

        fh = open('final.csv', 'w')
        fh.write(fileHeader)

        for iRegion in range(0, len(operationalRegions)):
            # Loop through all regionss
            curOperationalRegion = operationalRegions[iRegion]
            curSeverity = curOperationalRegion[0]
            curUpperDeltaV = curOperationalRegion[1]
            curUpperHeight = curOperationalRegion[2]
            curLowerDeltaV = curOperationalRegion[3]
            curLowerHeight = curOperationalRegion[4]

            # Create permutations of all possible height/speed combinations
            for iUpper in range(0, len(curUpperDeltaV)):
                for iLower in range(0, len(curLowerDeltaV)):
                    if curUpperHeight[iUpper] > curLowerHeight[iLower]:
                        severity.append(curSeverity)
                        upperHeight.append(curUpperHeight[iUpper])
                        upperDeltaV.append(curUpperDeltaV[iUpper])
                        lowerHeight.append(curLowerHeight[iLower])
                        lowerDeltaV.append(curLowerDeltaV[iLower])
                        result.append(RESULT_PENDING)
                        simulationFilename.append('')
                        areaTop.append(0)
                        batteryCapacity.append(0)
                        batteryWeight.append(0)

                        toProcess += 1

                        fh.write(str(curSeverity) + "; ")
                        fh.write(str(curUpperHeight[iUpper]) + "; ")
                        fh.write(str(curUpperDeltaV[iUpper]) + "; ")
                        fh.write(str(curLowerHeight[iLower]) + "; ")
                        fh.write(str(curLowerDeltaV[iLower]) + "; ")
                        fh.write(str(RESULT_PENDING) + "; ")
                        fh.write("nope.txt; ")
                        fh.write("0.0; ")
                        fh.write("0.0; ")
                        fh.write("0.0\n")

        fh.close()

    print(' >', toProcess, 'permutations left to process')
    ascentMaps = {}

    for iPermutation in range(0, len(severity)):
        # Check if the current solution wasn't already processed
        if result[iPermutation] != RESULT_PENDING:
            print('... skipping')
            continue

        # Check if an ascent map should be constructed
        lowerGuide = None

        if severity[iPermutation] in ascentMaps:
            # Already exists
            lowerGuide = ascentMaps[severity[iPermutation]]
        else:
            # Ascent guide does not exist yet
            print(TrackCommon.StringHeader("Determining ascent guides", 60))
            print(TrackCommon.StringPad(" > Severity:", severity[iPermutation], 2, 6))

            ascentGuide = TrackClimbOptimize.GenerateAscentMaps(axisHeight,
                axisDeltaV, settings.W, settings.S, settings.inclination,
                settings.lookupCl, settings.lookupCd, atmosphere, settings.qInfMin,
                settings.qInfMax, settings.alphaMin, settings.alphaMax, PReqMin, PReq,
                settings.latitude, settings.longitude, settings.reynoldsLength,
                severity=severity[iPermutation])
            
            ascentDeltaV = ascentGuide['pathMinDeltaV'] + climbLowerBoundOffset * \
                (ascentGuide['pathMaxDeltaV'] - ascentGuide['pathMinDeltaV'])

            lowerGuide = TrackLookup.Lookup1D(ascentGuide['pathHeight'],
                                              ascentDeltaV)

            ascentMaps[severity[iPermutation]] = lowerGuide

        settings.lowerBound = lowerGuide
        settings.upperBound = None

        # Start processing data
        print(TrackCommon.StringHeader("Simulating response", 60))
        print(TrackCommon.StringPad("severity:    ", severity[iPermutation], 3, 6))
        print(TrackCommon.StringPad("upper height:", upperHeight[iPermutation] / 1e3, 2, 6), "km")
        print(TrackCommon.StringPad("upper deltaV:", upperDeltaV[iPermutation], 2, 6), "m/s")
        print(TrackCommon.StringPad("lower height:", lowerHeight[iPermutation] / 1e3, 2, 6), "km")
        print(TrackCommon.StringPad("lower deltaV:", lowerDeltaV[iPermutation], 2, 6), "m/s")

        try:
            # Attempt to simulate track. It will throw an exception if
            # somehow the solution cannot be simulated or will not
            # converge
            stitched = TrackStitching.StitchTracks(
                upperHeight[iPermutation], upperDeltaV[iPermutation],
                lowerHeight[iPermutation], lowerDeltaV[iPermutation], 10,
                lowerDeltaV[iPermutation], upperDeltaV[iPermutation], PReqMin,
                PReq, dt, settings, severity[iPermutation], plotResult=False,
                saveResult=True)

            result[iPermutation] = RESULT_DONE
            simulationFilename[iPermutation] = stitched

            # Very backwards, I know. But I don't have the time nor the energy
            # and certainly not the motivation, at the moment, to do this
            # properly. And I have an SSD. Muhahaha
            loader = TrackStorage.DataStorage()
            loader.load(stitched)

            time = loader.getVariable('timeTotal').getValues()
            height = loader.getVariable('heightTotal').getValues()
            alpha = loader.getVariable('alphaTotal').getValues()
            gamma = loader.getVariable('gammaTotal').getValues()
            power = loader.getVariable('powerTotal').getValues()
            latitude = loader.getVariable('latitude').getValues()
            longitude = loader.getVariable('longitude').getValues()
            PRequired = loader.getVariable('PReqMax').getValues()

            # Determine area and capacity
            area, capacity = TrackStitching.DetermineArea(time, height, alpha,
                gamma, power, latitude, longitude, np.asarray([2.0]), PRequired,
                powerEstimator, settings, plotResults=False)

            batWeight, _ = AircraftBatteries.AircraftBatterySizing(
                capacity[0] / 3600, settings.batteryDepthOfDischarge,
                settings.batterySafetyFactor)

            areaTop[iPermutation] = area[0]
            batteryCapacity[iPermutation] = capacity[0]
            batteryWeight[iPermutation] = batWeight
        except Exception as ex:
            print(TrackCommon.StringHeader("ERROR: " + str(ex), 60, '!', '!', '!'))
            result[iPermutation] = RESULT_ERROR

        # Very inefficient. But the following means I can ctrl+c out of the
        # script at any time I want without issues with having to rerun
        # everything. Neat isn't it... Hacky and neat.
        print(TrackCommon.StringHeader("WRITING FILE, DO NOT CTRL+C", 60, '!', '!', '!'))
        fh = open(filename, 'w')

        fh.write(fileHeader)
        
        for i in range(0, len(severity)):
            fh.write(str(severity[i]) + "; ")
            fh.write(str(upperHeight[i]) + "; ")
            fh.write(str(upperDeltaV[i]) + "; ")
            fh.write(str(lowerHeight[i]) + "; ")
            fh.write(str(lowerDeltaV[i]) + "; ")
            fh.write(str(result[i]) + "; ")
            fh.write(simulationFilename[i] + "; ")
            fh.write(str(areaTop[i]) + "; ")
            fh.write(str(batteryCapacity[i]) + "; ")
            fh.write(str(batteryWeight[i]) + "\n")

        fh.close()

    # Save all final results to a csv file
    fh.close()

def PostProcess(filenames):
    settings = TrackSettings.Settings()
    atmosphere = Atmosphere.Atmosphere()

    for filename in filenames:
        file = TrackStorage.DataStorage()
        file.load(filename)

        powerEstimator = TrackPower.TrackPower(settings, atmosphere)

        time = file.getVariable('timeTotal').getValues()
        height = file.getVariable('heightTotal').getValues()
        alpha = file.getVariable('alphaTotal').getValues()
        gamma = file.getVariable('gammaTotal').getValues()
        power = file.getVariable('powerTotal').getValues()
        latitude = file.getVariable('latitude').getValues()
        longitude = file.getVariable('longitude').getValues()
        PRequired = file.getVariable('PReqMax').getValues()

        area, capacity = TrackStitching.DetermineArea(time, height, alpha,
            gamma, power, latitude, longitude, np.asarray([2.0]), PRequired,
            powerEstimator, settings, plotResults=False)

        area = area[0]
        capacity = capacity[0]

        print(' * For file:', filename)
        print(TrackCommon.StringPad(" > Atop:     ", area, 2, 6) + " m2")
        print(TrackCommon.StringPad(" > Capacity: ", capacity / 1e6, 2, 6) + " MJ")

        batWeight, batVolume = AircraftBatteries.AircraftBatterySizing(
            capacity / 3600, settings.batteryDepthOfDischarge,
            settings.batterySafetyFactor)

        solarPanelWeight = 2.0 * area * settings.specificWeightPanels

        print(TrackCommon.StringPad(" > mBattery: ", batWeight, 2, 6) + " kg")
        print(TrackCommon.StringPad(" > mSolar:   ", solarPanelWeight, 2, 6) + " kg")
        
def PostBatteryAnalysis(filename):
    file = TrackStorage.DataStorage()
    file.load(filename)
    numArea = 500
    numCapacities = 500
    
    time = file.getVariable('timeTotal').getValues()
    height = file.getVariable('heightTotal').getValues()
    alpha = file.getVariable('alphaTotal').getValues()
    gamma = file.getVariable('gammaTotal').getValues()
    power = file.getVariable('powerTotal').getValues()
    latitude = file.getVariable('latitude').getValues()
    longitude = file.getVariable('longitude').getValues()
    
    atmosphere = Atmosphere.Atmosphere()
    settings = TrackSettings.Settings()
    powerEstimator = TrackPower.TrackPower(settings, atmosphere)
    
    numSubUpdate = 10
    numUpdate = 50
    estimator = TimeEstimator.TimeEstimator(numArea)
    
    A = np.linspace(5, 40, numArea)
    Cmin = np.zeros(A.shape)
    Cmax = np.zeros(A.shape)
    
    # Precalculate constants related to efficiencies
    effDir = np.zeros(time.shape)
    effInd = np.zeros(time.shape)
    effUni = np.zeros(time.shape)
    
    for i in range(0, len(time)):
        effDir[i], effInd[i], effUni[i] = powerEstimator.getPowerEfficiency(
            height[i], latitude, longitude, alpha[i] / 180.0 * np.pi, gamma[i])
    
    estimator.startTiming()
    
    for iArea in range(0, len(A)):
        estimator.startIteration(iArea)
        
        isExcessArray = np.empty([len(time)], dtype=np.bool)
        isExcessArray.fill(False)
        
        curRequired = power / settings.efficiencyPropellers + settings.PConsumption
        curAvailable = 2 * A[iArea] * settings.fluxVenus * (effDir + effInd + 2.0 * effUni)

        # Turn the boolean list into a set of ranges
        isExcess = curAvailable[0] > curRequired[0]
        iOld = 0
        
        rangesExcess = []
        rangesLack = []

        for i in range(1, len(time)):
            if curAvailable[i] > curRequired[i]:
                if not isExcess:
                    # Current available is larger than required but was 
                    # previously tracking a non-excess region
                    rangesLack.append([iOld, i])
                    isExcess = True
                    iOld = i
            else:
                if isExcess:
                    # Currently lacking power but was tracking a charging
                    # region
                    rangesExcess.append([iOld, i])
                    isExcess = False
                    iOld = i
            
        # Add the last region
        if isExcess:
            rangesExcess.append([iOld, len(time)])
        else:
            rangesLack.append([iOld, len(time)])
            
        # Perform the intergral calculations
        upIntRequired = 0
        upIntAvailable = 0
        lowIntRequired = 0
        lowIntAvailable = 0
        
        for region in rangesExcess:
            upIntRequired += scp_int.simps(curRequired[region[0]:region[1]], time[region[0]:region[1]])
            upIntAvailable += scp_int.simps(curAvailable[region[0]:region[1]], time[region[0]:region[1]])
        
        for region in rangesLack:
            lowIntRequired += scp_int.simps(curRequired[region[0]:region[1]], time[region[0]:region[1]])
            lowIntAvailable += scp_int.simps(curAvailable[region[0]:region[1]], time[region[0]:region[1]])
            
        Cmax[iArea] = - settings.efficiencyCharging / settings.efficiencyPower * \
            upIntRequired + settings.efficiencyCharging * upIntAvailable
        Cmin[iArea] = 1 / (settings.efficiencyPower * settings.efficiencyCharging) * \
            lowIntRequired - lowIntAvailable / settings.efficiencyCharging
            
        estimator.finishedIteration(iArea)
        
        if (iArea + 1) % numSubUpdate == 0:
            print('.', end='')
        
        if (iArea + 1) % numUpdate == 0:
            print(' spent:', estimator.getTotalElapsed(), 
                  ', remaining:', estimator.getEstimatedRemaining())
            
    totalWeight = np.zeros([numCapacities, len(A)])
    capacities = np.linspace(np.min(Cmin), np.max(Cmax), numCapacities)
    
    for iArea in range(0, len(A)):
        for iVert in range(0, numCapacities):
            totalWeight[iVert, iArea], _ = AircraftBatteries.AircraftBatterySizing(
                capacities[iVert] / 3600, settings.batteryDepthOfDischarge, 
                settings.batterySafetyFactor)
            totalWeight[iVert, iArea] += A[iArea] * 2.0 * settings.specificWeightPanels * \
                settings.efficiencyPackingFactor

    # Put all excessive capacities to the maximum value
    maxWeight = np.max(totalWeight)
    
    for iArea in range(0, len(A)):
        for iVert in range(0, numCapacities):
            if capacities[iVert] < Cmin[iArea] or capacities[iVert] > Cmax[iArea]:
                totalWeight[iVert, iArea] = maxWeight
            
    boundsData, boundsLegend = TrackCommon.ImageAxes(0, 1, 0, 1)
    fig = plt.figure()
    ax = fig.add_axes(boundsData)
    axLegend = fig.add_axes(boundsLegend)
    TrackCommon.PlotImage(fig, ax, axLegend, A, r'$A_\mathrm{top} \; [m^2]$',
        capacities / 1e6, r'$C \; [MJ]$', totalWeight, r'$m_\mathrm{power} \; [kg]$')
    
    ax.plot(A, Cmax / 1e6, 'r', label=r'$C_\mathrm{max}$')
    ax.plot(A, Cmin / 1e6, 'g', label=r'$C_\mathrm{min}$')
    ax.legend()
    ax.grid(True)
            
#Main()
PostBatteryAnalysis('stitched_66076.9at5.0to50000.0at-22.2.dat')