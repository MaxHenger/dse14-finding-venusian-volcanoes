#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 15:53:48 2016

@author: MaxHenger
"""

import Atmosphere

import AircraftBatteries
import TrackCommon
import TrackClimbOptimize
import TrackAcceleratingOptimize
import TrackDiveOptimize
import TrackBounds
import TrackSettings
import TrackPower
import TrackContour
import TrackStorage
import TrackAngleOfAttack
import TrackLookup

import TimeEstimator

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scp_int

from mpl_toolkits.mplot3d import Axes3D

class MasterFile:
    def __init__(self):
        self.severities = []
        self.severityFiles = []
        self.simulationFiles = []

    def add(self, severity, severityFile='', simulationFile=''):
        self.severities.append(severity)
        self.severityFiles.append(severityFile)
        self.simulationFiles.append(simulationFile)
        
    def getNumSeverities(self):
        return len(self.severities)
    
    def getSeverityData(self, index):
        return self.severities[index], self.severityFiles[index], \
            self.simulationFiles[index]

    def save(self, filename):
        fh = open(filename, 'w')
        fh.write('severity; setup file; simulation file\n')
        
        for i in range(0, len(self.severities)):
            fh.write(str(self.severities[i]) + "; ")
            fh.write(self.severityFiles[i] + "; ")
            fh.write(self.simulationFiles[i] + "\n")
            
        fh.close()
        
    def load(self, filename):
        fh = open(filename, 'r')
        line = fh.readline() # discard header
        line = fh.readline()
        
        while len(line) != 0 and (len(line) != 1 or line[0] != '\n'):
            seperated = line.split(';')
            
            if len(seperated) != 3:
                raise ValueError("Expected three values per line")
            
            self.severities.append(float(seperated[0].strip()))
            self.severityFiles.append(seperated[1].strip())
            self.simulationFiles.append(seperated[2].strip())
            
            line = fh.readline()
            
        fh.close()
        
class SimulationFile:
    STATUS_PENDING = 0
    STATUS_DONECLIMB = 1
    STATUS_DONEDIVE = 2
    STATUS_DONE = 3
    STATUS_ERROR = 4
    
    def __init__(self, subfolder=''):
        if len(subfolder) != 0:
            if subfolder[-1] != '/':
                self.subfolder = subfolder + '/'
            else:
                self.subfolder = subfolder
        else:
            self.subfolder='./'
            
        self.lowerHeight = []
        self.lowerSpeedCruise = []
        self.lowerSpeedClimb = []
        self.upperHeight = []
        self.upperSpeedCruise = []
        self.lowerTimeCruise = []
        self.upperTimeCruise = []
        self.status = []
        self.dataDive = []
        self.dataAccPostDive = []
        self.dataLowerCruise = []
        self.dataAccPreClimb = []
        self.dataClimb = []
        self.dataAccPostClimb = []
        self.dataUpperCruise = []
        self.SPlanform = []
        self.SControl = []
        self.mSolarCell = []
        self.CBattery = []
        self.mBattery = []
        self.fOverlap = []
        self.mTotal = []

    def add(self, lowerHeight, lowerSpeedCruise, lowerSpeedClimb,  upperHeight, 
            upperSpeedCruise, lowerTimeCruise=0, upperTimeCruise=0, status=STATUS_PENDING,
            dataDive='', dataAccPostDive='', dataLowerCruise='', dataAccPreClimb='',
            dataClimb='', dataAccPostClimb='', dataUpperCruise='', SPlanform=0,
            SControl=0, mSolarCell=0, CBattery=0, mBattery=0, fOverlap=-1, 
            mTotal=0):
        self.lowerHeight.append(lowerHeight)
        self.lowerSpeedCruise.append(lowerSpeedCruise)
        self.lowerSpeedClimb.append(lowerSpeedClimb)
        self.upperHeight.append(upperHeight)
        self.upperSpeedCruise.append(upperSpeedCruise)
        self.lowerTimeCruise.append(lowerTimeCruise)
        self.upperTimeCruise.append(upperTimeCruise)
        self.status.append(status)
        self.dataDive.append(dataDive)
        self.dataAccPostDive.append(dataAccPostDive)
        self.dataLowerCruise.append(dataLowerCruise)
        self.dataAccPreClimb.append(dataAccPreClimb)
        self.dataClimb.append(dataClimb)
        self.dataAccPostClimb.append(dataAccPostClimb)
        self.dataUpperCruise.append(dataUpperCruise)
        self.SPlanform.append(SPlanform)
        self.SControl.append(SControl)
        self.mSolarCell.append(mSolarCell)
        self.CBattery.append(CBattery)
        self.mBattery.append(mBattery)
        self.fOverlap.append(fOverlap)
        self.mTotal.append(mTotal)
        
    def getNumEntries(self):
        return len(self.lowerHeight)
        
    def save(self, filename):
        fh = open(filename, 'w')
        
        fh.write('h_low [m]; v_cr_low [m/s]; v_cl_low [m/s]; h_up [m]; v_cr_up [m/s]; ' +
            't_cr_low [s]; t_cr_up [s]; status; data_dive; data_acc_post_dive; ' +
            'data_cr_low; data_acc_pre_climb; data_climb; data_acc_post_climb; ' +
            'data_cr_up; S_plan [m2]; S_cont [m2]; m_sc [kg]; C_bat [J]; m_bat [kg]; ' + 
            'f_overlap; m_tot [kg]\n')
        
        for i in range(0, len(self.lowerHeight)):
            fh.write(str(self.lowerHeight[i]) + '; ')
            fh.write(str(self.lowerSpeedCruise[i]) + '; ')
            fh.write(str(self.lowerSpeedClimb[i]) + '; ')
            fh.write(str(self.upperHeight[i]) + '; ')
            fh.write(str(self.upperSpeedCruise[i]) + '; ')
            fh.write(str(self.lowerTimeCruise[i]) + '; ')
            fh.write(str(self.upperTimeCruise[i]) + '; ')
            fh.write(str(self.status[i]) + '; ')
            fh.write(str(self.dataDive[i]) + '; ')
            fh.write(str(self.dataAccPostDive[i]) + '; ')
            fh.write(str(self.dataLowerCruise[i]) + '; ')
            fh.write(str(self.dataAccPreClimb[i]) + '; ')
            fh.write(str(self.dataClimb[i]) + '; ')
            fh.write(str(self.dataAccPostClimb[i]) + '; ')
            fh.write(str(self.dataUpperCruise[i]) + '; ')
            fh.write(str(self.SPlanform[i]) + '; ')
            fh.write(str(self.SControl[i]) + '; ')
            fh.write(str(self.mSolarCell[i]) + '; ')
            fh.write(str(self.CBattery[i]) + '; ')
            fh.write(str(self.mBattery[i]) + '; ')
            fh.write(str(self.fOverlap[i]) + '; ')
            fh.write(str(self.mTotal[i]) + '\n')
            
        fh.close()
            
    def load(self, filename):
        fh = open(filename, 'r')
        
        line = fh.readline() # header
        line = fh.readline()
        
        while len(line) != 0 and (len(line) != 1 or line[0] != '\n'):
            seperated = line.split(';')
            
            if len(seperated) != 22:
                raise ValueError("Expected 22 values per line")
                
            self.lowerHeight.append(float(seperated[0].strip()))
            self.lowerSpeedCruise.append(float(seperated[1].strip()))
            self.lowerSpeedClimb.append(float(seperated[2].strip()))
            self.upperHeight.append(float(seperated[3].strip()))
            self.upperSpeedCruise.append(float(seperated[4].strip()))
            self.lowerTimeCruise.append(float(seperated[5].strip()))
            self.upperTimeCruise.append(float(seperated[6].strip()))
            self.status.append(int(seperated[7].strip()))
            self.dataDive.append(seperated[8].strip())
            self.dataAccPostDive.append(seperated[9].strip())
            self.dataLowerCruise.append(seperated[10].strip())
            self.dataAccPreClimb.append(seperated[11].strip())
            self.dataClimb.append(seperated[12].strip())
            self.dataAccPostClimb.append(seperated[13].strip())
            self.dataUpperCruise.append(seperated[14].strip())
            self.SPlanform.append(float(seperated[15].strip()))
            self.SControl.append(float(seperated[16].strip()))
            self.mSolarCell.append(float(seperated[17].strip()))
            self.CBattery.append(float(seperated[18].strip()))
            self.mBattery.append(float(seperated[19].strip()))
            self.fOverlap.append(float(seperated[20].strip()))
            self.mTotal.append(float(seperated[21].strip()))
            
            line = fh.readline()
            
        fh.close()
        
def DetermineOperatingPoints(axisHeight, axisDeltaV, axisSeverity, latitude,
                             longitude, numUpper, numLower, filename):
    # Settings for the script 
    severityBoundaries = 0.1 # percent to stay away from minimum and maximum severities
    numSeverities = 3
    deltaVBoundaries = 0.1
    deltaVPerOp = 15
    minHeightOp = 32e3
    maxHeightOp = 52e3
    deltaHeightPerOp = 5e3
    numTestHeight = 100
    
    # Create the managing classes used
    settings = TrackSettings.Settings()
    atmosphere = Atmosphere.Atmosphere()
    powerEstimator = TrackPower.TrackPower(settings, atmosphere)
    timeEstimator = TimeEstimator.TimeEstimator(len(axisSeverity))
    
    lookupdCldAlpha = settings.lookupCl.getDerivative()
    lookupdCddAlpha = settings.lookupCd.getDerivative()
    
    STopTotal = settings.S + settings.SControl
    
    # Settings for this specific script    
    print(TrackCommon.StringHeader("Determine operational points", 60))
    
    # Loop through all provided severities
    timeEstimator.startTiming()
    
    severityMin = None
    severityMax = None
    
    for iSeverity in range(len(axisSeverity)):
        timeEstimator.startIteration(iSeverity)
        
        # Generate a cruise map and check if sustained cruise is possible
        cruiseMaps = TrackBounds.GenerateCruiseMaps(axisHeight, axisDeltaV, 
            settings.latitude, settings.longitude, axisSeverity[iSeverity],
            settings.W, settings.S, settings.inclination, settings.lookupCl,
            settings.lookupCd, atmosphere, settings.reynoldsLength)
        
        cruiseContour = TrackContour.Contour()
        cruiseContour.combineData(axisDeltaV, axisHeight, [cruiseMaps['qInf'],
            cruiseMaps['alpha'], cruiseMaps['PReq'], cruiseMaps['vInfOvera']],
            [settings.qInfMin, settings.alphaMin, 0, 0], [settings.qInfMax,
            settings.alphaMax, settings.PPropeller, settings.speedOfSoundRatio])
        
        # Generate structures required for storage
        sustainedCruiseHeight = []
        sustainedCruiseVHor = []
        sustainedCruisePProp = []
        sustainedCruisePAvailable = []
        sustainedCruiseAlpha = []
        
        # Find the contour that contains the main cruising region
        mainContour = None
        
        for iContour in range(cruiseContour.getNumContours()):
            curContour = cruiseContour.getContour(iContour)
            
            if len(curContour.getVerticalRanges(0)) != 0:
                mainContour = curContour
                break
            
        if mainContour == None:
            if severityMin != None:
                severityMax = axisSeverity[iSeverity - 1]
                break
        else:
            # Find if a cruise speed exists where sustained cruise is possible
            allowedSustainedCruise = False
            
            minContourHeight = mainContour.getMinY()
            maxContourHeight = mainContour.getMaxY()
            
            iMinContourHeight = TrackCommon.Find1DBisectionAscending(axisHeight, minContourHeight)
            iMaxContourHeight = TrackCommon.Find1DBisectionAscending(axisHeight, maxContourHeight)
            
            iFirstValidCruise = None
            iLastValidCruise = None
            
            for iCruiseHeight in range(iMinContourHeight + 1, iMaxContourHeight):
                # For the current height determine the horizontal velocity that
                # is required to stay at the same subsolar point, then figure
                # out if it is a point within the cruise contours
                curVHor = settings.omegaVenus * (settings.RVenus + axisHeight[iCruiseHeight])
                cruiseRanges = mainContour.getVerticalRanges(curVHor)
                
                isContained = False
                
                for iRange in range(len(cruiseRanges)):
                    if cruiseRanges[iRange][0] <= axisHeight[iCruiseHeight] and \
                            cruiseRanges[iRange][1] >= axisHeight[iCruiseHeight]:
                        isContained = True
                        break
                    
                if not isContained:
                    continue
                
                # Determine the required thrust and angle of attack
                curDensity = TrackCommon.AdjustSeverity(atmosphere.density(
                    axisHeight[iCruiseHeight], settings.latitude, settings.longitude),
                    axisSeverity[iSeverity])
                    
                curVZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(
                    axisHeight[iCruiseHeight], settings.latitude, settings.longitude),
                    axisSeverity[iSeverity])
                    
                curVInf = curVZonal + curVHor
                curQInf = 0.5 * curDensity * curVInf**2.0
                
                curRe = TrackCommon.AdjustSeverity(atmosphere.reynoldsNumber(
                    axisHeight[iCruiseHeight], settings.latitude, settings.longitude,
                    curVInf, settings.reynoldsLength), axisSeverity[iSeverity])
            
                alpha, thrust, valid = TrackAngleOfAttack.AngleOfAttackThrustSteady(
                    settings.W, settings.S, curQInf, settings.inclination, 
                    settings.lookupCl, lookupdCldAlpha, settings.lookupCd, 
                    lookupdCddAlpha, curRe)
                
                if not valid:
                    continue
                
                PRequired = thrust * curVInf / settings.efficiencyPropellers + settings.PConsumption
                
                # Determine the available power
                effDir, effInd, effUni = powerEstimator.getPowerEfficiency(
                    axisHeight[iCruiseHeight], settings.latitude, 
                    settings.longitude, alpha / 180.0 * np.pi, 0)
                
                PAvailable = settings.fluxVenus * STopTotal * (effDir + effInd + 2.0 * effUni)
                
                if PAvailable * settings.efficiencyPower > PRequired:
                    allowedSustainedCruise = True
                    sustainedCruiseHeight.append(axisHeight[iCruiseHeight])
                    sustainedCruiseVHor.append(curVHor)
                    sustainedCruisePProp.append(thrust * curVInf)
                    sustainedCruisePAvailable.append(PAvailable)
                    sustainedCruiseAlpha.append(alpha)

            if allowedSustainedCruise and severityMin is None:
                severityMin = axisSeverity[iSeverity]
        
        timeEstimator.finishedIteration(iSeverity)
        print('spent:', timeEstimator.getTotalElapsed(),
              ', remaining:', timeEstimator.getEstimatedRemaining())
        
    print(' > min severity:', round(severityMin, 3))
    print(' > max severity:', round(severityMax, 3))
    
    axisSeverity = np.linspace(severityMin + severityBoundaries * 
        (severityMax - severityMin) / 2, severityMax - severityBoundaries *
        (severityMax - severityMin) / 2, numSeverities)
    severityFilenames = []
    
    for iSeverity in range(numSeverities):
        # For the current severity determine the operating points and the
        # ascent maps
        data = TrackStorage.DataStorage()
        
        ascentMaps = TrackClimbOptimize.GenerateAscentMaps(axisHeight,
            axisDeltaV, settings.W, settings.S, settings.inclination, 
            settings.lookupCl, settings.lookupCd, atmosphere, settings.qInfMin, 
            settings.qInfMax, settings.alphaMin, settings.alphaMax, 0, 
            settings.PPropeller, settings.latitude, settings.longitude, 
            settings.reynoldsLength, axisSeverity[iSeverity], storeResults=False)
        cruiseMaps = TrackBounds.GenerateCruiseMaps(axisHeight, axisDeltaV,
            settings.latitude, settings.longitude, axisSeverity[iSeverity],
            settings.W, settings.S, settings.inclination, settings.lookupCl, 
            settings.lookupCd, atmosphere, settings.reynoldsLength)
        
        ascentPathDelta = ascentMaps['pathMaxDeltaV'] - ascentMaps['pathMinDeltaV']
        ascentPathMean = (ascentMaps['pathMaxDeltaV'] + ascentMaps['pathMinDeltaV']) / 2.0
        
        
        cruiseContour = TrackContour.Contour()
        cruiseContour.combineData(axisDeltaV, axisHeight, [cruiseMaps['qInf'],
            cruiseMaps['alpha'], cruiseMaps['PReq'], cruiseMaps['vInfOvera']],
            [settings.qInfMin, settings.alphaMin, 0, 0], [settings.qInfMax,
            settings.alphaMax, settings.PPropeller, settings.speedOfSoundRatio])
        
        # Find the proper cruising contour
        mainCruiseContour = None
        
        for i in range(cruiseContour.getNumContours()):
            curContour = cruiseContour.getContour(i)
            
            if len(curContour.getVerticalRanges(0)):
                mainCruiseContour = curContour
                break
            
        rangeDeltaV = mainCruiseContour.getMaxX()
        numDeltaV = int(rangeDeltaV * (1.0 - deltaVBoundaries) / deltaVPerOp) + 1
        upperOperationalDeltaV = None
        upperOperationalHeight = []
        upperOperationalExcess = []
        
        if numDeltaV == 1:
            upperOperationalDeltaV = [rangeDeltaV / 2]
        else:
            offsetDeltaV = deltaVBoundaries * rangeDeltaV / 2
            upperOperationalDeltaV = np.linspace(offsetDeltaV, rangeDeltaV - 
                offsetDeltaV, numDeltaV)
        
        # For each of the operational deltaV values find the corresponding 
        # height at which the ratio of available power to required power is
        # maximized
        for iOp in range(numDeltaV):
            curDeltaV = upperOperationalDeltaV[iOp]
            heightRanges = mainCruiseContour.getVerticalRanges(curDeltaV)
            
            bestExcess = None
            bestHeight = 0
            
            for iRange in range(len(heightRanges)):
                curRange = heightRanges[iRange]

                testHeight = np.linspace(curRange[0], curRange[1], numTestHeight)
                vZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(
                    testHeight, settings.latitude, settings.longitude), 
                    axisSeverity[iSeverity])
                density = TrackCommon.AdjustSeverity(atmosphere.density(
                    testHeight, settings.latitude, settings.longitude),
                    axisSeverity[iSeverity])
                vInf = vZonal + curDeltaV
                qInf = 0.5 * density * np.power(vInf, 2)
                
                for iHeight in range(numTestHeight):
                    Re = TrackCommon.AdjustSeverity(atmosphere.reynoldsNumber(
                        testHeight[iHeight], settings.latitude, settings.longitude,
                        vInf[iHeight], settings.reynoldsLength), 
                        axisSeverity[iSeverity])
                
                    alpha, thrust, valid = TrackAngleOfAttack.AngleOfAttackThrustSteady(
                        settings.W, settings.S, qInf[iHeight], settings.inclination,
                        settings.lookupCl, lookupdCldAlpha, settings.lookupCd,
                        lookupdCddAlpha, Re)
                    
                    effDir, effInd, effUni = powerEstimator.getPowerEfficiency(
                        testHeight[iHeight], settings.latitude, settings.longitude,
                        alpha / 180.0 * np.pi, 0)
                    
                    PExcess = STopTotal * settings.fluxVenus * (effDir + effInd +
                        2.0 * effUni) * settings.efficiencyPower - (thrust * vInf[iHeight] / 
                        settings.efficiencyPropellers + settings.PConsumption)
                    
                    if bestExcess is None or PExcess > bestExcess:
                        bestExcess = PExcess
                        bestHeight = testHeight[iHeight]

            upperOperationalHeight.append(bestHeight)
            upperOperationalExcess.append(bestExcess)
            
        # Find the lower segment operational heights
        lookupPreAscent = TrackLookup.Lookup1D(ascentMaps['axisHeight'],
            ascentPathMean + 0.5 * ascentPathDelta)
        
        lowerMinHeight = minHeightOp
        lowerMaxHeight = maxHeightOp
        
        heightRanges = mainCruiseContour.getVerticalRanges(0)
        
        if mainCruiseContour.getMinY() > lowerMinHeight:
            lowerMinHeight = mainCruiseContour.getMinY()
            
        if len(heightRanges) != 0 and heightRanges[-1][1] < lowerMaxHeight:
            lowerMaxHeight = heightRanges[-1][1]

        numLowerOp = int((lowerMaxHeight - lowerMinHeight) / deltaHeightPerOp) + 1

        lowerOperationalHeight = list(np.linspace(lowerMinHeight, lowerMaxHeight, numLowerOp))
        lowerOperationalDeltaV = []
        lowerOperationalMaxDeltaV = []
        lowerToRemove = []
        
        for iHeight in range(numLowerOp):
            curDeltaV = lookupPreAscent(lowerOperationalHeight[iHeight])
            lowerOperationalDeltaV.append(curDeltaV)
            
            if curDeltaV > 0:
                lowerToRemove.append(iHeight)
                
            # Retrieve maximum cruise speed at bottom
            speedRanges = mainCruiseContour.getHorizontalRanges(lowerOperationalHeight[iHeight])
            lowerOperationalMaxDeltaV.append(speedRanges[-1][-1])
                
        for iRemove in range(len(lowerToRemove)):
            del lowerOperationalHeight[lowerToRemove[iRemove] - iRemove]
            del lowerOperationalDeltaV[lowerToRemove[iRemove] - iRemove]
            del lowerOperationalMaxDeltaV[lowerToRemove[iRemove] - iRemove]
            
        # Find the height and deltaV values at which continuous cruise is 
        # possible
        iMinContourHeight = TrackCommon.Find1DBisectionAscending(axisHeight,
            mainCruiseContour.getMinX())
        iMaxContourHeight = TrackCommon.Find1DBisectionAscending(axisHeight,
            mainCruiseContour.getMaxX())
        
        constantCruiseHeight = []
        constantCruiseDeltaV = []
        
        for iHeight in range(iMinContourHeight + 1, iMaxContourHeight):
            curVHor = settings.omegaVenus * (settings.RVenus + axisHeight[iCruiseHeight])
            cruiseRanges = mainCruiseContour.getVerticalRanges(curVHor)
            
            # Check if the current height is contained by the ranges
            # provided by the current vHor
            isContained = False
            for iRange in range(len(cruiseRanges)):
                if cruiseRanges[iRange][0] <= axisHeight[iHeight] and \
                        cruiseRanges[iRange][1] >= axisHeight[iHeight]:
                    isContained = True
                    break
                
            if not isContained:
                continue
            
            # Calculate all atmospheric and flight properties at the currently
            # considered height and deltaV
            curDensity = TrackCommon.AdjustSeverity(atmosphere.density(
                axisHeight[iHeight], settings.latitude, settings.longitude),
                axisSeverity[iSeverity])
            curVZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(
                axisHeight[iHeight], settings.latitude, settings.longitude),
                axisSeverity[iSeverity])
            curVInf = curVZonal + curVHor
            curQInf = 0.5 * density * curVInf**2.0
            curRe = TrackCommon.AdjustSeverity(atmosphere.reynoldsNumber(
                axisHeight[iHeight], settings.latitude, settings.longitude,
                curVInf, settings.reynoldsLength), axisSeverity[iSeverity])
            
            alpha, thrust, valid = TrackAngleOfAttack.AngleOfAttackThrustSteady(
                settings.W, settings.S, curQInf, settings.inclination,
                settings.lookupCl, lookupdCldAlpha, settings.lookupCd,
                lookupdCddAlpha, curRe)
            
            if not valid:
                continue
            
            # Calculate if there is any excess power, if so then this cruising
            # point is considered valid
            effDir, effInd, effUni = powerEstimator.getPowerEfficiency(
                axisHeight[iHeight], settings.latitude, settings.longitude,
                alpha / 180.0 * np.pi, 0)
            
            if settings.fluxVenus * STopTotal * (effDir + effInd + 2.0 * effUni) * \
                    settings.efficiencyPower > thrust * curVInf / \
                    settings.efficiencyPropellers + settings.PConsumption:
                constantCruiseHeight.append(axisHeight[iHeight])
                constantCruiseDeltaV.append(curVHor)
            
        # Store all retrieved data in a storage file
        data.addVariable('axisHeight', axisHeight)
        data.addVariable('axisDeltaV', axisDeltaV)
        data.addVariable('severity', axisSeverity[iSeverity])
        data.addVariable('ascentGuideHeight', ascentMaps['axisHeight'])
        data.addVariable('ascentGuideDeltaV', ascentPathMean)
        data.addVariable('upperOperationalHeight', upperOperationalHeight)
        data.addVariable('upperOperationalDeltaV', upperOperationalDeltaV)
        data.addVariable('lowerOperationalHeight', lowerOperationalHeight)
        data.addVariable('lowerOperationalDeltaV', lowerOperationalDeltaV)
        data.addVariable('lowerOperationalMaxDeltaV', lowerOperationalMaxDeltaV)
        data.addVariable('constantCruiseHeight', constantCruiseHeight)
        data.addVariable('constantCruiseDeltaV', constantCruiseDeltaV)
        data.addVariable('ascentVVer', ascentMaps['vVer'])
        data.addVariable('ascentPReq', ascentMaps['PReq'])
        data.addVariable('ascentAlpha', ascentMaps['alpha'])
        data.addVariable('ascentQInf', ascentMaps['qInf'])
        data.addVariable('cruiseQInf', cruiseMaps['qInf'])
        data.addVariable('cruiseVInf', cruiseMaps['vInf'])
        data.addVariable('cruiseAlpha', cruiseMaps['alpha'])
        data.addVariable('cruisePReq', cruiseMaps['PReq'])
        data.addVariable('latitude', settings.latitude)
        data.addVariable('longitude', settings.longitude)
        data.addVariable('W', settings.W)
        data.addVariable('S', settings.S)
        data.addVariable('PPropeller', settings.PPropeller)
        data.addVariable('speedOfSoundRatio', settings.speedOfSoundRatio)
        data.addVariable('lookupClAxisX', settings.lookupCl.getPoints()[0])
        data.addVariable('lookupClAxisY', settings.lookupCl.getPoints()[1])
        data.addVariable('lookupCl', settings.lookupCl.getPoints()[2])
        data.addVariable('lookupCdAxisX', settings.lookupCd.getPoints()[0])
        data.addVariable('lookupCdAxisY', settings.lookupCd.getPoints()[1])
        data.addVariable('lookupCd', settings.lookupCd.getPoints()[2])
        data.addVariable('reynoldsLength', settings.reynoldsLength)
        data.addVariable('inclination', settings.inclination)
        curFilename = 'severity_' + str(round(axisSeverity[iSeverity], 3))
        data.save(curFilename)
        severityFilenames.append(curFilename)
        
    masterFile = MasterFile()

    for i in range(len(severityFilenames)):
        masterFile.add(axisSeverity[i], severityFilenames[i])
    
    masterFile.save(filename)
    
def SaveSimulationSection(time, height, vHor, vVer, vInf, alpha, gamma, power, filename):
    data = TrackStorage.DataStorage()
    data.addVariable('time', time)
    data.addVariable('height', height)
    data.addVariable('vHor', vHor)
    data.addVariable('vVer', vVer)
    data.addVariable('vInf', vInf)
    data.addVariable('alpha', alpha)
    data.addVariable('gamma', gamma)
    data.addVariable('power', power)
    data.save(filename)
    
def LoadSimulationSection(filename):
    data = TrackStorage.DataStorage()
    data.load(filename)
    time = data.getVariable('time').getValues()
    height = data.getVariable('height').getValues()
    vHor = data.getVariable('vHor').getValues()
    vVer = data.getVariable('vVer').getValues()
    vInf = data.getVariable('vInf').getValues()
    alpha = data.getVariable('alpha').getValues()
    gamma = data.getVariable('gamma').getValues()
    power = data.getVariable('power').getValues()
    
    return time, height, vHor, vVer, vInf, alpha, gamma, power
    
def AppendSimulationSection(time, height, vHor, vVer, vInf, alpha, gamma, power, filename):
    nTime, nHeight, nVHor, nVVer, nVInf, nAlpha, nGamma, nPower = \
        LoadSimulationSection(filename)
        
    time.extend(nTime)
    height.extend(nHeight)
    vHor.extend(nVHor)
    vVer.extend(nVVer)
    vInf.extend(nVInf)
    alpha.extend(nAlpha)
    gamma.extend(nGamma)
    power.extend(nPower)
    
def CalculatePowerIntegrals(time, height, alpha, gamma, power, latitude,
                            longitude, settings, A, estimator):
    PRequired = power / settings.efficiencyPropellers + settings.PConsumption
    PAvailable = np.zeros(time.shape)
    
    for i in range(len(time)):
        effDir, effInd, _ = estimator.getPowerEfficiency(height[i],
            latitude, longitude, alpha[i] / 180.0 * np.pi, gamma[i])
        
        PAvailable[i] = A * (effDir + effInd) * settings.fluxVenus

    charging = []
    discharging = []

    isCharging = PAvailable[0] * settings.efficiencyPower > PRequired[0]
    iOld = 0
    
    for i in range(1, len(time)):
        if PAvailable[i] * settings.efficiencyPower > PRequired[i]:
            if not isCharging:
                # Go into charging mode
                discharging.append([iOld, i + 1])
                iOld = i
                isCharging = True
        elif isCharging:
            # Go into discharging mode
            charging.append([iOld, i + 1])
            iOld = i
            isCharging = False
            
    if isCharging:
        charging.append([iOld, len(time)])
    else:
        discharging.append([iOld, len(time)])
            
    # Calculate integrals associated with the charging and discharging regions
    intChargeRequired = 0
    intChargeAvailable = 0
    intDischargeRequired = 0
    intDischargeAvailable = 0
    
    for curRange in charging:
        intChargeRequired += scp_int.simps(PRequired[curRange[0]:curRange[1]],
                                           time[curRange[0]:curRange[1]])
        intChargeAvailable += scp_int.simps(PAvailable[curRange[0]:curRange[1]],
                                            time[curRange[0]:curRange[1]])
        
    for curRange in discharging:
        intDischargeRequired += scp_int.simps(PRequired[curRange[0]:curRange[1]],
                                              time[curRange[0]:curRange[1]])
        intDischargeAvailable += scp_int.simps(PAvailable[curRange[0]:curRange[1]],
                                               time[curRange[0]:curRange[1]])
    
    return intChargeRequired, intChargeAvailable, \
        intDischargeRequired, intDischargeAvailable
        
def GenerateTracks(mainFile, severity, dt, tol=1e-3, forceRedive=False, forceRepower=False):
    # Settings for this script
    boundsDeltaV = 0.1
    cruiseDeltaVStep = 5
    subFolder = './trackResults/v2/'
    atmosphere = Atmosphere.Atmosphere()
    settings = TrackSettings.Settings()
    powerEstimator = TrackPower.TrackPower(settings, atmosphere)
    
    # Load the main file and its contents
    masterFile = MasterFile()
    masterFile.load(mainFile)
    
    iMasterFile = 0
    boundsFilename = ''
    simulationFilename = ''
    
    for i in range(masterFile.getNumSeverities()):
        curSeverity, curBoundsFilename, curSimulationFilename = masterFile.getSeverityData(i)
        
        if abs(curSeverity - severity) < tol:
            iMasterFile = i
            severity = curSeverity
            boundsFilename = curBoundsFilename
            simulationFilename = curSimulationFilename
            
    if len(boundsFilename) == 0:
        raise ValueError("Failed to find specified severity in master file")
        
    # Load the bounds file
    data = TrackStorage.DataStorage()
    data.load(boundsFilename)
    
    upperHeight = data.getVariable('upperOperationalHeight').getValues()
    upperDeltaV = data.getVariable('upperOperationalDeltaV').getValues()
    lowerHeight = data.getVariable('lowerOperationalHeight').getValues()
    lowerDeltaV = data.getVariable('lowerOperationalDeltaV').getValues()
    lowerMaxDeltaV = data.getVariable('lowerOperationalMaxDeltaV').getValues()
    ascentGuideHeight = data.getVariable('ascentGuideHeight').getValues()
    ascentGuideDeltaV = data.getVariable('ascentGuideDeltaV').getValues()
    
    lookupCl = TrackLookup.Lookup2D(data.getVariable('lookupClAxisX').getValues(),
                                    data.getVariable('lookupClAxisY').getValues(),
                                    data.getVariable('lookupCl').getValues())
    lookupCd = TrackLookup.Lookup2D(data.getVariable('lookupCdAxisX').getValues(),
                                    data.getVariable('lookupCdAxisY').getValues(),
                                    data.getVariable('lookupCd').getValues())
    lookupdCldAlpha = lookupCl.getDerivative()
    lookupdCddAlpha = lookupCd.getDerivative()
    
    ascentGuide = TrackLookup.Lookup1D(ascentGuideHeight, ascentGuideDeltaV)
    
    latitude = data.getVariable('latitude').getValues()
    longitude = data.getVariable('longitude').getValues()
    W = data.getVariable('W').getValues()
    S = data.getVariable('S').getValues()
    PPropeller = data.getVariable('PPropeller').getValues()
    speedOfSoundRatio = data.getVariable('speedOfSoundRatio').getValues()
    reynoldsLength = data.getVariable('reynoldsLength').getValues()
    inclination = data.getVariable('inclination').getValues()
    STop = (S + settings.SControl) * settings.efficiencyPackingFactor
    
    # Check if a simulation file is already present or if one should be 
    # generated
    simulationFile = SimulationFile(subFolder)
    
    if len(simulationFilename) == 0:
        # Need to generate a simulation file
        print(TrackCommon.StringHeader("Generating simulation file", 60))
        
        for iLower in range(len(lowerHeight)):
            # Determine the speeds at which one can fly in the lower track
            deltaLowerSpeed = lowerMaxDeltaV[iLower] - lowerDeltaV[iLower]
            numLowerSpeeds = int(deltaLowerSpeed * (1 - boundsDeltaV) / cruiseDeltaVStep) + 1
            lowerSpeeds = None
            
            if numLowerSpeeds == 1:
                lowerSpeeds = [(lowerMaxDeltaV[iLower] + lowerDeltaV[iLower]) / 2.0]
            else:
                lowerSpeeds = np.linspace(lowerDeltaV[iLower] + 0.5 * boundsDeltaV * deltaLowerSpeed,
                    lowerMaxDeltaV[iLower] - 0.5 * boundsDeltaV * deltaLowerSpeed, numLowerSpeeds)
            
            for iLowerSpeeds in range(numLowerSpeeds):
                for iUpper in range(len(upperHeight)):
                    simulationFile.add(lowerHeight[iLower], lowerSpeeds[iLowerSpeeds],
                        lowerDeltaV[iLower], upperHeight[iUpper], upperDeltaV[iUpper])
                    
        simulationFilename = 'sim_' + str(round(severity, 3)) + '.csv'
        masterFile.simulationFiles[iMasterFile] = simulationFilename
        
        print(TrackCommon.StringHeader("Saving file, do not press <Ctrl+C>", 60, 'S', 'S', 'S'))
        simulationFile.save(simulationFilename)
        masterFile.save(mainFile)
        print(TrackCommon.StringHeader("Done saving, safe to quit", 60, 'S', 'S', 'S'))
    else:
        print(TrackCommon.StringHeader("Loading simulation file", 60))
        simulationFile.load(simulationFilename)
        
    # With the simulation file loaded first check for all simulations that
    # still require the climb track to be generated
    for iSimulation in range(simulationFile.getNumEntries()):
        if simulationFile.status[iSimulation] == simulationFile.STATUS_PENDING:
            # Perform climbing simulation
            deltaHeight = simulationFile.upperHeight[iSimulation] - \
                simulationFile.lowerHeight[iSimulation]
            climbFilename = None
            climbStatus = None
            
            try:
                time, height, vHor, vVer, vInf, alpha, gamma = \
                    TrackClimbOptimize.OptimizeClimb(simulationFile.lowerHeight[iSimulation],
                        simulationFile.upperHeight[iSimulation], simulationFile.lowerHeight[iSimulation] -
                        0.1 * deltaHeight, simulationFile.lowerSpeedClimb[iSimulation], 0, longitude,
                        latitude, W, S, simulationFile.upperSpeedCruise[iSimulation], 0, PPropeller,
                        speedOfSoundRatio, inclination, dt, lookupCl, lookupCd, reynoldsLength,
                        severity=severity, lookupBoundLowerVInf=ascentGuide, plotResults=False,
                        storeResults=False)
                    
                climbFilename = 'climb_' + str(round(severity, 3)) + '_' + \
                    str(round(simulationFile.lowerHeight[iSimulation], 3)) + '_to_' + \
                    str(round(simulationFile.upperHeight[iSimulation], 3)) + '.dat'
                SaveSimulationSection(time, height, vHor, vVer, vInf, alpha, 
                    gamma, np.repeat(PPropeller, len(time)), subFolder + climbFilename)
                climbStatus = simulationFile.STATUS_DONECLIMB
            except Exception as ex:
                print(TrackCommon.StringHeader("EXCEPTION: " + str(ex), 60, '!', '!', '!'))
                
                climbFilename = 'ERROR'
                climbStatus = simulationFile.STATUS_ERROR
            
            for iOther in range(simulationFile.getNumEntries()):
                if abs(simulationFile.lowerHeight[iSimulation] - simulationFile.lowerHeight[iOther]) < tol and \
                        abs(simulationFile.upperHeight[iSimulation] - simulationFile.upperHeight[iOther]) < tol and \
                        abs(simulationFile.lowerSpeedClimb[iSimulation] - simulationFile.lowerSpeedClimb[iOther]) < tol:
                    simulationFile.dataClimb[iOther] = subFolder + climbFilename
                    simulationFile.status[iOther] = climbStatus

            print(TrackCommon.StringHeader("Saving file, do not press <Ctrl+C>", 60, 'S', 'S', 'S'))
            simulationFile.save(simulationFilename)
            print(TrackCommon.StringHeader("Done saving, safe to quit", 60, 'S', 'S', 'S'))
    
    # Second round: check which files need to have the accelerating, cruising,
    # and diving phases simulated
    for iSimulation in range(simulationFile.getNumEntries()):
        if simulationFile.status[iSimulation] == simulationFile.STATUS_DONECLIMB or forceRedive:
            # Load the climbing file and figure out the final height and 
            # horizontal speed
            dataClimb = TrackStorage.DataStorage()
            dataClimb.load(simulationFile.dataClimb[iSimulation])
            heightClimb = dataClimb.getVariable('height').getValues()
            vHorClimb = dataClimb.getVariable('vHor').getValues()
            alphaClimb = dataClimb.getVariable('alpha').getValues()
            
            # Accelerate to cruising speed
            accPostClimbFilename = None
            hasError = False
            
            avgVHor = 0
            avgPower = 0
            
            try:
                time, vHor, alpha, power, avgVHor, avgPower = \
                    TrackAcceleratingOptimize.OptimizeAccelerating(heightClimb[-1],
                        vHorClimb[-1], PPropeller, alphaClimb[-1], longitude,
                        latitude, W, S, simulationFile.upperSpeedCruise[iSimulation],
                        inclination, dt, 0, PPropeller, speedOfSoundRatio, lookupCl,
                        lookupCd, reynoldsLength, severity=severity, plotResults=False,
                        storeResults=False)
                    
                time.append(time[-1] + dt)
                vHor.append(avgVHor)
                alpha.append(alpha[-1])
                power.append(avgPower)
                
                accPostClimbFilename = "accPostClimb_" + str(round(severity, 3)) + '_' + \
                    str(round(simulationFile.upperHeight[iSimulation], 3)) + '_' + \
                    str(round(vHorClimb[-1], 3)) + 'to' + \
                    str(round(simulationFile.upperSpeedCruise[iSimulation], 3)) + '.dat'
                    
                vZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(
                    heightClimb[-1], latitude, longitude), severity)
                    
                SaveSimulationSection(time, np.repeat(heightClimb[-1], len(time)),
                    vHor, np.repeat(0.0, len(time)), vHor + vZonal, alpha,
                    np.repeat(0.0, len(time)), power, subFolder + accPostClimbFilename)
                
            except Exception as ex:
                print(TrackCommon.StringHeader("EXCEPTION: " + str(ex), 60, '!', '!', '!'))
                
                accPostClimbFilename = 'ERROR'
                hasError = True
            
            simulationFile.dataAccPostClimb[iSimulation] = subFolder + accPostClimbFilename

            if hasError:
                simulationFile.status[iSimulation] = simulationFile.STATUS_ERROR
                
                print(TrackCommon.StringHeader("Saving file, do not press <Ctrl+C>", 60, 'S', 'S', 'S'))
                simulationFile.save(simulationFilename)
                print(TrackCommon.StringHeader("Done saving, safe to quit", 60, 'S', 'S', 'S'))
                
                continue
            
            # Perform dive
            diveFilename = None
            heightDive = 0
            vHorDive = 0
            alphaDive = 0
            
            try:
                time, height, vHor, vVer, vInf, alpha, gamma = \
                    TrackDiveOptimize.OptimizeDive(heightClimb[-1], heightClimb[0],
                        avgVHor, 0, longitude, latitude, W, S, simulationFile.lowerSpeedCruise[iSimulation],
                        0, speedOfSoundRatio, dt, lookupCl, lookupCd, reynoldsLength, severity,
                        plotResults=False, storeResults=False)
                    
                diveFilename = 'dive_' + str(round(severity, 3)) + '_' + \
                    str(round(simulationFile.upperHeight[iSimulation], 3)) + 'at' + \
                    str(round(simulationFile.upperSpeedCruise[iSimulation], 3)) + '_' + \
                    str(round(simulationFile.lowerHeight[iSimulation], 3)) + 'at' + \
                    str(round(simulationFile.lowerSpeedCruise[iSimulation], 3)) + '.dat'
                    
                SaveSimulationSection(time, height, vHor, vVer, vInf, alpha,
                    gamma, 0.0, subFolder + diveFilename)
                
                heightDive = height[-1]
                vHorDive = vHor[-1]
                alphaDive = alpha[-1]
                
            except Exception as ex:
                print(TrackCommon.StringHeader("EXCEPTION: " + str(ex), 60, '!', '!', '!'))
                
                diveFilename = 'ERROR'
                hasError = True
                
            simulationFile.dataDive[iSimulation] = subFolder + diveFilename

            if hasError:
                simulationFile.status[iSimulation] = simulationFile.STATUS_ERROR

                print(TrackCommon.StringHeader("Saving file, do not press <Ctrl+C>", 60, 'S', 'S', 'S'))
                simulationFile.save(simulationFilename)
                print(TrackCommon.StringHeader("Done saving, safe to quit", 60, 'S', 'S', 'S'))
                
                continue
            
            # Perform post-dive acceleration
            accPostDiveFilename = None
            accAlpha = 0
            
            try:
                time, vHor, alpha, power, avgVHor, avgPower = \
                    TrackAcceleratingOptimize.OptimizeAccelerating(heightDive,
                        vHorDive, 0, alphaDive, longitude, latitude, W, S,
                        simulationFile.lowerSpeedCruise[iSimulation],
                        inclination, dt, 0, PPropeller, speedOfSoundRatio,
                        lookupCl, lookupCd, reynoldsLength, severity=severity,
                        plotResults=False, storeResults=False)
                    
                time.append(time[-1] + dt)
                vHor.append(avgVHor)
                alpha.append(alpha[-1])
                power.append(avgPower)
                
                accAlpha = alpha[-1]
                
                accPostDiveFilename = "accPostDive_" + str(round(severity, 3)) + '_' + \
                    str(round(simulationFile.lowerHeight[iSimulation], 3)) + '_' + \
                    str(round(vHorDive, 3)) + 'to' + \
                    str(round(simulationFile.lowerSpeedCruise[iSimulation], 3)) + '.dat'
                    
                vZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(
                    heightDive, latitude, longitude), severity)
                
                SaveSimulationSection(time, np.repeat(heightDive, len(time)),
                    vHor, np.repeat(0.0, len(time)), vHor + vZonal, alpha,
                    np.repeat(0.0, len(time)), power, subFolder + accPostDiveFilename)
                
            except Exception as ex:
                print(TrackCommon.StringHeader("EXCEPTION: " + str(ex), 60, '!', '!', '!'))
                
                accPostDiveFilename = 'ERROR'
                hasError = True
                
            simulationFile.dataAccPostDive[iSimulation] = subFolder + accPostDiveFilename
            
            if hasError:
                simulationFile.status[iSimulation] = simulationFile.STATUS_ERROR

                print(TrackCommon.StringHeader("Saving file, do not press <Ctrl+C>", 60, 'S', 'S', 'S'))
                simulationFile.save(simulationFilename)
                print(TrackCommon.StringHeader("Done saving, safe to quit", 60, 'S', 'S', 'S'))
                
                continue
            
            # Perform post-cruise acceleration to climb speed
            accPreClimbFilename = None
            
            try:
                oldVHor = avgVHor
                time, vHor, alpha, power, avgVHor, avgPower = \
                    TrackAcceleratingOptimize.OptimizeAccelerating(heightDive,
                        avgVHor, avgPower, accAlpha, longitude, latitude, W, S,
                        simulationFile.lowerSpeedClimb[iSimulation],
                        inclination, dt, 0, PPropeller, speedOfSoundRatio,
                        lookupCl, lookupCd, reynoldsLength, severity,
                        plotResults=False, storeResults=False)
                    
                time.append(time[-1] + dt)
                vHor.append(avgVHor)
                alpha.append(alpha[-1])
                power.append(avgPower)
                
                accPreClimbFilename = "accPreDive_" + str(round(severity, 3)) + "_" + \
                    str(round(simulationFile.lowerHeight[iSimulation], 3)) + '_' + \
                    str(round(oldVHor, 3)) + 'to' + \
                    str(round(simulationFile.lowerSpeedClimb[iSimulation], 3)) + '.dat'
                    
                vZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(
                    heightDive, latitude, longitude), severity)
                
                SaveSimulationSection(time, np.repeat(heightDive, len(time)),
                    vHor, np.repeat(0.0, len(time)), vHor + vZonal, alpha,
                    np.repeat(0.0, len(time)), power, subFolder + accPreClimbFilename)
                
            except Exception as ex:
                print(TrackCommon.StringHeader("EXCEPTION: " + str(ex), 60, '!', '!', '!'))
                
                accPreClimbFilename = 'ERROR'
                hasError = True
                
            simulationFile.dataAccPreClimb[iSimulation] = subFolder + accPreClimbFilename

            if hasError:
                simulationFile.status[iSimulation] = simulationFile.STATUS_ERROR

                print(TrackCommon.StringHeader("Saving file, do not press <Ctrl+C>", 60, 'S', 'S', 'S'))
                simulationFile.save(simulationFilename)
                print(TrackCommon.StringHeader("Done saving, safe to quit", 60, 'S', 'S', 'S'))
                
                continue
            
            # Successfully processed all accelerating and diving phases, 
            # proceed to the power phase with the current file
            simulationFile.status[iSimulation] = simulationFile.STATUS_DONEDIVE
            
            print(TrackCommon.StringHeader("Saving file, do not press <Ctrl+C>", 60, 'S', 'S', 'S'))
            simulationFile.save(simulationFilename)
            print(TrackCommon.StringHeader("Done saving, safe to quit", 60, 'S', 'S', 'S'))
            
    # Third round, check which files have complete tracks but still require the
    # power-related parameters and the cruising times to be estimated
    CDesign = AircraftBatteries.AircraftCapacity(settings.mBattery, 
        settings.batteryDepthOfDischarge, settings.batterySafetyFactor)
    
    # Generate test capacities (defined in mass, then converted to )
    mTestBatteries = np.linspace(50, 250, 10)
    CTest = np.zeros(mTestBatteries.shape)
    
    for i in range(len(mTestBatteries)):
        CTest[i] = AircraftBatteries.AircraftCapacity(mTestBatteries[i],
            settings.batteryDepthOfDischarge, settings.batterySafetyFactor)
        
    mTestBatteries = np.append(mTestBatteries, settings.mBattery)
    CTest = np.append(CTest, CDesign)
    numLowerTime = 30
    
    totalMBattery = []
    totalTimeLower = []
    totalTimeUpper = []
    totalHeightLower = []
    totalHeightUpper = []
    totalFOverlap = []
    
    for iSimulation in range(simulationFile.getNumEntries()):
        if simulationFile.status[iSimulation] == simulationFile.STATUS_DONEDIVE or \
                (forceRepower and simulationFile.status[iSimulation] == simulationFile.STATUS_DONE):
            timeLower = []
            timeUpper = []
            heightLower = []
            heightUpper = []
            fOverlap = []
            valid = []

            # Load all simulation data
            timeDive, heightDive, vHorDive, vVerDive, vInfDive, alphaDive, \
                gammaDive, powerDive = LoadSimulationSection(simulationFile.dataDive[iSimulation])
                
            timePostDive, heightPostDive, vHorPostDive, vVerPostDive, \
                vInfPostDive, alphaPostDive, gammaPostDive, powerPostDive = \
                LoadSimulationSection(simulationFile.dataAccPostDive[iSimulation])
            
            timePreClimb, heightPreClimb, vHorPreClimb, vVerPreClimb, \
                vInfPreClimb, alphaPreClimb, gammaPreClimb, powerPreClimb, = \
                LoadSimulationSection(simulationFile.dataAccPreClimb[iSimulation])
                
            timeClimb, heightClimb, vHorClimb, vVerClimb, vInfClimb, \
                alphaClimb, gammaClimb, powerClimb = \
                LoadSimulationSection(simulationFile.dataClimb[iSimulation])
                
            timePostClimb, heightPostClimb, vHorPostClimb, vVerPostClimb, \
                vInfPostClimb, alphaPostClimb, gammaPostClimb, powerPostClimb = \
                LoadSimulationSection(simulationFile.dataAccPostClimb[iSimulation])
                
            # Calculate all integrals associated with these phases
            intDiveChR, intDiveChA, intDiveDisR, intDiveDisA = CalculatePowerIntegrals(
                timeDive, heightDive, alphaDive, gammaDive, np.repeat(0.0, len(timeDive)), latitude,
                longitude, settings, STop, powerEstimator)
            
            intPostDiveChR, intPostDiveChA, intPostDiveDisR, intPostDiveDisA = \
                CalculatePowerIntegrals(timePostDive, heightPostDive, alphaPostDive, 
                gammaPostDive, powerPostDive, latitude, longitude, settings,
                STop, powerEstimator)
                
            intPreClimbChR, intPreClimbChA, intPreClimbDisR, intPreClimbDisA = \
                CalculatePowerIntegrals(timePreClimb, heightPreClimb, alphaPreClimb,
                gammaPreClimb, powerPreClimb, latitude, longitude, settings,
                STop, powerEstimator)
                
            intClimbChR, intClimbChA, intClimbDisR, intClimbDisA = \
                CalculatePowerIntegrals(timeClimb, heightClimb, alphaClimb, gammaClimb,
                powerClimb, latitude, longitude, settings, STop, powerEstimator)
                
            intPostClimbChR, intPostClimbChA, intPostClimbDisR, intPostClimbDisA = \
                CalculatePowerIntegrals(timePostClimb, heightPostClimb, alphaPostClimb, 
                gammaPostClimb, powerPostClimb, latitude, longitude, settings,
                STop, powerEstimator)
                
            totalChR = intDiveChR + intPostDiveChR + intPreClimbChR + \
                intClimbChR + intPostClimbChR
            totalChA = intDiveChA + intPostDiveChA + intPreClimbChA + \
                intClimbChA + intPostClimbChA
            totalDisR = intDiveDisR + intPostDiveDisR + intPreClimbDisR + \
                intClimbDisR + intPostClimbDisR
            totalDisA = intDiveDisA + intPostDiveDisA + intPreClimbDisA + \
                intClimbDisA + intPostClimbDisA
                
            # Calculate cruising quantities
            # - lower cruising
            lowerZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(
                simulationFile.lowerHeight[iSimulation], latitude, longitude), severity)
            lowerDensity = TrackCommon.AdjustSeverity(atmosphere.density(
                simulationFile.lowerHeight[iSimulation], latitude, longitude), severity)
            lowerVInf = lowerZonal + simulationFile.lowerSpeedCruise[iSimulation]
            lowerQInf = 0.5 * lowerDensity * lowerVInf**2.0
            lowerRe = TrackCommon.AdjustSeverity(atmosphere.reynoldsNumber(
                simulationFile.lowerHeight[iSimulation], latitude, longitude,
                lowerVInf, reynoldsLength), severity)
            lowerAlpha, lowerThrust, lowerValid = TrackAngleOfAttack.AngleOfAttackThrustSteady(
                W, S, lowerQInf, inclination, lookupCl, lookupdCldAlpha, lookupCd, 
                lookupdCddAlpha, lowerRe)
            
            if not lowerValid:
                print('not valid lower!')
                continue
            
            lowerEfficiencies = powerEstimator.getPowerEfficiency(
                simulationFile.lowerHeight[iSimulation], latitude, 
                longitude, lowerAlpha / 180.0 * np.pi, 0)
            lowerPower = lowerThrust * lowerVInf
            lowerR = lowerPower / settings.efficiencyPropellers + settings.PConsumption
            lowerA = settings.fluxVenus * STop * (lowerEfficiencies[0] + lowerEfficiencies[1])
            
            # - upper cruising
            upperZonal = TrackCommon.AdjustSeverity(atmosphere.velocityZonal(
                simulationFile.upperHeight[iSimulation], latitude, longitude), severity)
            upperDensity = TrackCommon.AdjustSeverity(atmosphere.density(
                simulationFile.upperHeight[iSimulation], latitude, longitude), severity)
            upperVInf = upperZonal + simulationFile.upperSpeedCruise[iSimulation]
            upperQInf = 0.5 * upperDensity * upperVInf**2.0
            upperRe = TrackCommon.AdjustSeverity(atmosphere.reynoldsNumber(
                simulationFile.upperHeight[iSimulation], latitude, longitude,
                upperVInf, reynoldsLength), severity)
            upperAlpha, upperThrust, upperValid = TrackAngleOfAttack.AngleOfAttackThrustSteady(
                W, S, upperQInf, inclination, lookupCl, lookupdCldAlpha, lookupCd, 
                lookupdCddAlpha, upperRe)
            
            if not upperValid:
                print('not valid upper!')
                continue
            
            upperEfficiencies = powerEstimator.getPowerEfficiency(
                simulationFile.upperHeight[iSimulation], latitude, longitude,
                upperAlpha / 180.0 * np.pi, 0)
            upperPower = upperThrust * upperVInf
            upperR = upperPower / settings.efficiencyPropellers + settings.PConsumption
            upperA = settings.fluxVenus * STop * (upperEfficiencies[0] + upperEfficiencies[1])
   
            Cmin = (1 / (settings.efficiencyDischarging * settings.efficiencyPower) * \
                    (totalDisR - totalDisA * settings.efficiencyPower))
            
            # Loop through all battery sizes
            for iBattery in range(len(CTest)):
                # Check if the aircraft is able to fly at the given track by 
                # calculating the time the aircraft can spend in the lower track
                newTimeLower = []
                newTimeUpper = []
                newHeightUpper = []
                newHeightLower = []
                newFOverlap = []

                if CTest[iBattery] < Cmin:
                    continue
                
                tLowerMax = (settings.efficiencyDischarging * settings.efficiencyPower *
                    CTest[iBattery] - (intDiveDisR + intPostDiveDisR + intPreClimbDisR +
                    intClimbDisR + intPostClimbDisR) + (intDiveDisA + intPostDiveDisA + 
                    intPreClimbDisA + intClimbDisA + intPostClimbDisA) *
                    settings.efficiencyPower) / (lowerR - 
                    settings.efficiencyPower * lowerA)
                    
                if tLowerMax < 0:
                    # Invalid solution
                    continue
                elif tLowerMax > 1800:
                    # ADD TEMPERATURE ESTIMATION CODE HERE
                    tLowerMax = 1800
                
                # Create a range of lower cruising times
                tLower = np.linspace(0, tLowerMax, numLowerTime)
                
                for iLowerTime in range(numLowerTime):
                    curTLower = tLower[iLowerTime]

                    # Calculate the upper time required to provide the battery
                    # with enough energy to perform the lower phases
                    tUpper = ( 1 / (settings.efficiencyDischarging * settings.efficiencyPower) *
                        ( intDiveDisR + intPostDiveDisR + intPreClimbDisR + intClimbDisR +
                        intPostClimbDisR - settings.efficiencyPower * ( intDiveDisA + intPostDiveDisA +
                        intPreClimbDisA + intClimbDisA + intPostClimbDisA ) + curTLower * 
                        (lowerR - settings.efficiencyPower * lowerA)) - settings.efficiencyPower *
                        settings.efficiencyCharging * ( intDiveChA + intPostDiveChA + 
                        intPreClimbChA + intClimbChA + intPostClimbChA - (1 / settings.efficiencyPower) * 
                        (intDiveChR + intPostDiveChR + intPreClimbChR + intClimbChR + intPostClimbChR))) / \
                        (settings.efficiencyPower * settings.efficiencyCharging * (upperA - 
                        upperR / settings.efficiencyPower))
                    
                    # Calculate the minimum overlap length
                    intDiveGr = scp_int.simps(vHorDive / (settings.RVenus + heightDive) * 
                        settings.RVenus, timeDive)
                    intPostDiveGr = scp_int.simps(vHorPostDive / (settings.RVenus +
                        heightPostDive) * settings.RVenus, timePostDive)
                    intPreClimbGr = scp_int.simps(vHorPreClimb / (settings.RVenus +
                        heightPreClimb) * settings.RVenus, timePreClimb)
                    intClimbGr = scp_int.simps(vHorClimb / (settings.RVenus + 
                        heightClimb) * settings.RVenus, timeClimb)
                    intPostClimbGr = scp_int.simps(vHorPostClimb / (settings.RVenus +
                        heightPostClimb) * settings.RVenus, timePostClimb)
                    
                    d_bottom = intPostDiveGr + intPreClimbGr + curTLower * (
                        simulationFile.lowerSpeedCruise[iSimulation] * settings.RVenus /
                        (settings.RVenus + simulationFile.lowerHeight[iSimulation]))
                    
                    d_total = d_bottom + intClimbGr + intPostClimbGr + intDiveGr + \
                        tUpper * (simulationFile.upperSpeedCruise[iSimulation] *
                        settings.RVenus / (settings.RVenus + simulationFile.lowerHeight[iSimulation]))
                        
                    curFOverlap = d_total / d_bottom
                
                    newTimeLower.append(curTLower)
                    newTimeUpper.append(tUpper)
                    newHeightUpper.append(simulationFile.upperHeight[iSimulation])
                    newHeightLower.append(simulationFile.lowerHeight[iSimulation])
                    newFOverlap.append(curFOverlap)
                    
#                    print(' * Solution:')
#                    print(' > tLower: ', curTLower)
#                    print(' > tUpper: ', tUpper)
#                    print(' > hLower: ', simulationFile.lowerHeight[iSimulation])
#                    print(' > hUpper: ', simulationFile.upperHeight[iSimulation])
#                    print(' > overlap:', curFOverlap)
                    
                totalMBattery.append(mTestBatteries[iBattery])
                totalTimeLower.append(newTimeLower)
                totalTimeUpper.append(newTimeUpper)
                totalHeightLower.append(newHeightLower)
                totalHeightUpper.append(newHeightUpper)
                totalFOverlap.append(newFOverlap)
        
    if len(totalTimeLower) != 0:
        cmap = plt.get_cmap('gnuplot2')
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(np.asarray(totalHeightLower) / 1e3, np.asarray(totalHeightUpper) / 1e3, totalFOverlap,
                   totalMBattery)
        #ax.set_zlim(0, 10)
            
            
#DetermineOperatingPoints(np.linspace(20e3, 80e3, 121), np.linspace(-80, 80, 65),
#    np.linspace(-2.0, 2.0, 11), 0, 0, 5, 5, 'master.csv')
GenerateTracks('master.csv', 0.0, 0.35)
#GenerateTracks('master.csv', -1.44, 0.35)
#GenerateTracks('master.csv', 1.44, 0.35)