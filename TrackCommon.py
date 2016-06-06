# -*- coding: utf-8 -*-
"""
Created on Sun May 29 19:25:25 2016

@author: MaxHenger
"""

import numpy as np
import TrackLookup

import matplotlib as mpl
import matplotlib.pyplot as plt

def IsArray(var):
    return isinstance(var, list) or isinstance(var, np.ndarray)

def StringPad(name, value, decimals, length):
    valueString = str(round(value, decimals))

    dot = valueString.find('.', 0)

    if dot != -1:
        if len(valueString) - dot - 1 < decimals:
            valueString = valueString + ('0' * (decimals + 1 + dot - len(valueString)))

    if len(valueString) < length:
        valueString = ' ' * (length - len(valueString)) + valueString

    return name + valueString

def StringHeader(name, width, knot='+', hor='-', ver='|'):
    stringUpper = knot + ((width - 2) * hor) + knot
    stringCenter = ''

    if len(name) <= width - 4:
        remaining = int(width - 4 - len(name))

        if remaining % 2 == 0:
            stringCenter = ver + (' ' * int(remaining / 2 + 1)) + name + (' ' * int(remaining / 2 + 1)) + ver
        else:
            stringCenter = ver + (' ' * (int(remaining / 2) + 2)) + \
                name + (' ' * (int(remaining / 2) + 1)) + ver

    return stringUpper + '\n' + stringCenter + '\n' + stringUpper

def ImageAxes(xl, xr, yb, yt, colorbarSeperation=0.05, colorbarWidth=0.03, border=0.1):
    return [ # The plot's axes
        xl + border / 2,
        yb + border / 2,
        xr - xl - colorbarWidth - colorbarSeperation - border,
        yt - yb - border
    ], \
    [ # The colorbar's axes
        xr - colorbarWidth - colorbarSeperation / 2 - border / 2,
        yb + border / 2,
        colorbarWidth,
        yt - yb - border
    ]

def PlotImage(fig, axImage, axColorbar, xAxis, xLabel, yAxis, yLabel, data, dataLabel,
              cmap='gnuplot2', contours=None, forceNormMin=None, forceNormMax=None):
    # Determine extent of axes
    xMin = min(xAxis)
    xMax = max(xAxis)
    yMin = min(yAxis)
    yMax = max(yAxis)

    # Determine normalization
    vMin = 0
    vMax = 0

    if forceNormMin == None:
        vMin = np.min(data)
    else:
        vMin = forceNormMin

    if forceNormMax == None:
        vMax = np.max(data)
    else:
        vMax = forceNormMax

    norm = mpl.colors.Normalize(vMin, vMax)

    # Determine colormap
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)

    # Determine extent. shift by half a pixel leftward and downward if required.
    # This will assume all data is equally spaced
    extent = [xMin, xMax, yMin, yMax]
    axImage.imshow(data, extent=extent, norm=norm,
                   aspect='auto', origin='lower', cmap=cmap)
    axImage.set_xlabel(xLabel)
    axImage.set_ylabel(yLabel)
    axImage.grid(True)

    if contours != None:
        ct = axImage.contour(xAxis, yAxis, data, contours, colors='k')
        axImage.clabel(ct)

    # Add the colorbar
    cbb = mpl.colorbar.ColorbarBase(axColorbar, cmap=cmap, norm=norm)
    cbb.set_label(dataLabel, rotation=90, fontsize=14)

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

def AdjustSeverity(atmosphere, severity):
    if severity > 0.0:
        return atmosphere[1] + (atmosphere[2] - atmosphere[1]) * severity

    return atmosphere[1] - (atmosphere[1] - atmosphere[0]) * severity
    
def AdjustBiasMapIndividually(biasMap, amount, location, width, name):
    print('adjusting', name, 'bias by', round(amount, 3), 'at h', round(location, 1), 'km')

    if not biasMap.modifyCentered(amount, location, width):
        base = biasMap.getBaseBias()
        print('bias became negative, adjusting', name, 'base from',
            round(base, 3), 'to', round(base - amount, 3))

        biasMap.reset(base - amount)

    return biasMap.getBaseBias()

def AdjustBiasMapCommonly(biasMaps, amount, names):
    newBiases = [0] * len(biasMaps)

    for iMap in range(0, len(biasMaps)):
        curBias = biasMaps[i].getBaseBias()
        print('adjusting', names[i], 'bias from', round(curBias, 3),
            'to', round(curBias - amount, 3))

        newBiases[iMap] = curBias - amount
        biasMap.reset(newBiases[iMap])

    return newBiases

def Lerp(x1, y1, x2, y2, x):
    return y1 + (y2 - y1) / (x2 - x1) * (x - x1)

def __TestLoadAerodynamicData__():
    Cl, Cd = LoadAerodynamicData("./data/aerodynamicPerformance/Cl.csv",
                                 "./data/aerodynamicPerformance/Cd.csv")

    ClPoints = Cl.getPoints()
    CdPoints = Cd.getPoints()

    for i in range(0, len(ClPoints[0])):
        print('alpha =', round(ClPoints[0][i], 4),
              ', Cl =', round(ClPoints[1][i], 4),
              ', Cd =', round(CdPoints[1][i], 4))

#__TestLoadAerodynamicData__()
