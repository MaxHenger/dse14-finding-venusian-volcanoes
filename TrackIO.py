#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 00:21:47 2016

@author: MaxHenger
"""

import TrackLookup
import TrackStorage

import numpy as np

def LoadAerodynamicData(dataCl, dataCd, tol=1e-15):
    # Load data from file
    dataCl = np.genfromtxt(dataCl, delimiter=';')
    dataCd = np.genfromtxt(dataCd, delimiter=';')

    # Check for consistency of dimensions
    if dataCl.shape[1] < 2:
        raise ValueError("Expected at least two columns in Cl data file")

    if dataCd.shape[1] < 2:
        raise ValueError("Expected at least two columns in Cd data file")

    if dataCl.shape[0] != dataCd.shape[0]:
        raise ValueError("Expected same number of rows in Cl as in Cd file")

    # Check for consistency of angles of attack
    for iAlpha in range(0, dataCl.shape[0]):
        if abs(dataCl[iAlpha, 0] - dataCd[iAlpha, 0]) > tol:
            raise ValueError('ClAlpha =', dataCl[iAlpha, 0], '!=',
                             'CdAlpha =', dataCd[iAlpha, 0], 'at index', iAlpha)

    # Create the lookup tables and return them
    return TrackLookup.Lookup1D(dataCl[:, 0], dataCl[:, 1]), \
        TrackLookup.Lookup1D(dataCd[:, 0], dataCd[:, 1])

def _LoadAerodynamicReynoldsFile_(filename, tol=1e-15):
    # Open up the specified file
    fh = open(filename, 'r')

    line = fh.readline()

    # Read through the header and retrieve the positions at which the angles of
    # attack and the reynolds numbers are set
    iAlpha = []
    iReynolds = []
    reynolds = []

    headerEntries = line.split(',')

    for iEntry in range(0, len(headerEntries)):
        curEntry = headerEntries[iEntry].lower().strip()

        if curEntry.find('alpha') != -1:
            iAlpha.append(iEntry)
            continue

        indexRe = curEntry.find('re')

        if indexRe != -1:
            iReynolds.append(iEntry)
            curEntry = curEntry[indexRe + 2:]
            indexE = curEntry.find('e')

            if indexE == -1:
                raise ValueError("Found a column entry without an 'e' in the " +
                    "Reynolds number designation")

            preE = curEntry[:indexE]
            postE = curEntry[indexE + 1:]

            reynolds.append(float(preE) * 10 ** float(postE))

    iReynoldsOrdering = np.argsort(reynolds)
    reynolds = np.asarray(reynolds)[iReynoldsOrdering]

    if len(iAlpha) != len(iReynolds) or len(iAlpha) != len(reynolds):
        raise RuntimeError("Number of found alpha columns did not match the " +
            "number of found Reynolds columns")

    # Load data systematically. Keep track of the start and end where all
    # datapoints are valid (as a 2D matrix has to be constructed)
    iFirstValid = -1 # inclusive index
    iLastValid = -1 # exclusive index

    alpha = []
    values = []

    counter = 0

    line = fh.readline()

    while len(line) != 0 and (len(line) > 1 or line[-1] != '\n'):
        # Create a new list with angle-of-attack values
        values.append([])

        # Prepare the current line for reading. Check whether all values are
        # valid to determine when the final valid line of data is encountered
        seperated = line.split(',')
        allValid = True

        for iOrdered in iReynoldsOrdering:
            # Check the angle of attack and reynolds number for validity and
            # store them if they're both valid
            stringAlpha = seperated[iAlpha[iOrdered]].strip()
            stringValue = seperated[iReynolds[iOrdered]].strip()

            if len(stringValue) == 0 or len(stringAlpha) == 0:
                allValid = False
                values[-1].append(0.0)
                continue

            curAlpha = float(stringAlpha)

            if len(alpha) == counter:
                # First new alpha entry
                alpha.append(curAlpha)
            else:
                # New alpha entry, ensure it doens't differ by too much
                if abs(alpha[-1] - curAlpha) > tol:
                    raise RuntimeError("Encountered alpha differs too much from " +
                        "already stored alpha")

            values[-1].append(float(stringValue))

        # Check the validity of all datapoints
        if allValid:
            if iFirstValid == -1:
                iFirstValid = counter
        elif iFirstValid != -1 and iLastValid == -1:
            iLastValid = counter

        counter += 1
        line = fh.readline()

    fh.close()

    # Create the 2D lookup table
    return TrackLookup.Lookup2D(alpha[iFirstValid:iLastValid],
        reynolds[iFirstValid:iLastValid], values[iFirstValid:iLastValid])

def LoadAerodynamicReynoldsData(dataCl, dataCd, tol=1e-15):
    # Load data from file
    ClLookup = _LoadAerodynamicReynoldsFile_(dataCl, tol)
    CdLookup = _LoadAerodynamicReynoldsFile_(dataCd, tol)
    
    # Ensure the angles of attack and the reynolds numbers for both axes match
    axisAlphaCl, axisReCl, _ = ClLookup.getPoints()
    axisAlphaCd, axisReCd, _ = CdLookup.getPoints()
    
    if len(axisAlphaCl) != len(axisAlphaCd):
        raise ValueError("Expected the number of alpha values to match for the " +
            "Cl lookup table and the Cd lookup table")
        
    if len(axisReCl) != len(axisReCd):
        raise ValueError("Expected the number of Reynolds number values to " +
            "match for the Cl lookup table and the Cd lookup table")
        
    for i in range(0, len(axisAlphaCl)):
        if abs(axisAlphaCl[i] - axisAlphaCd[i]) > tol:
            raise ValueError("Expected the alpha arrays for the Cl lookup table " +
                "and the Cd lookup table to match")
            
    for i in range(0, len(axisReCl)):        
        if abs(axisReCl[i] - axisReCd[i]) > tol:
            raise ValueError("Expected the Reynolds number arrays for the Cl lookup " +
                "table and the Cd lookup table to match")
            
    return ClLookup, CdLookup

def LoadAscentGuides(data, safety=0.0):
    file = TrackStorage.DataStorage()
    file.load(data)

    height = file.getVariable('height').getValues()
    minimum = file.getVariable('minDeltaV').getValues()
    maximum = file.getVariable('maxDeltaV').getValues()

    lookupMinimum = TrackLookup.Lookup1D(height, minimum)
    lookupMaximum = TrackLookup.Lookup1D(height, maximum + (maximum - minimum) * safety)

    return lookupMinimum, lookupMaximum

def __TestLoadAerodynamicData__():
    Cl, Cd = LoadAerodynamicData("./data/aerodynamicPerformance/Cl.csv",
                                 "./data/aerodynamicPerformance/Cd.csv")

    ClPoints = Cl.getPoints()
    CdPoints = Cd.getPoints()

    for i in range(0, len(ClPoints[0])):
        print('alpha =', round(ClPoints[0][i], 4),
              ', Cl =', round(ClPoints[1][i], 4),
              ', Cd =', round(CdPoints[1][i], 4))

import matplotlib.pyplot as plt
def __TestLoadAerodynamicReynoldsData__():
    Cl, Cd = LoadAerodynamicReynoldsData('./data/aerodynamicPerformance/v2/ClCruiser.csv',
                                         './data/aerodynamicPerformance/v2/CdCruiser.csv')

    cmap = plt.get_cmap('YlOrRd')
    fig = plt.figure()
    
    # Plot Cl data
    axCl = fig.add_subplot(131)

    axisAlpha, axisRe, valuesCl = Cl.getPoints()

    for i in range(0, len(axisRe)):
        axCl.plot(axisAlpha, valuesCl[:, i], color=cmap(float(i + 6) / (len(axisRe) + 5)))
        
    # Plot Cd data
    axCd = fig.add_subplot(132)
    
    axisAlpha, axisRe, valuesCd = Cd.getPoints()
    
    for i in range(0, len(axisRe)):
        axCd.plot(axisAlpha, valuesCd[:, i], color=cmap(float(i + 6) / (len(axisRe) + 5)))
    
    # Plot cl/cd data
    valuesClOverCd = np.zeros(valuesCl.shape)
    
    for i in range(0, len(axisAlpha)):
        for j in range(0, len(axisRe)):
            valuesClOverCd[i, j] = valuesCl[i, j] / valuesCd[i, j]

    axClOverCd = fig.add_subplot(133)
    
    lines = []
    names = []
    
    for i in range(0, len(axisRe)):
        curLine, = axClOverCd.plot(axisAlpha, valuesClOverCd[:, i],
            color=cmap(float(i + 6) / (len(axisRe) + 5)))
        lines.append(curLine)
        names.append(r'$Re\;=\;' + str(round(axisRe[i] / 1e6, 3)) + r'$')
        
    axClOverCd.legend(lines, names, loc='best')

#__TestLoadAerodynamicData__()
#__TestLoadAerodynamicReynoldsData__()
