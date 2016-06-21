# -*- coding: utf-8 -*-
"""
Created on Sun May 29 19:25:25 2016

@author: MaxHenger
"""

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

# IsArray will return true if the variable passed into it is an array (that is:
# a numpy array or a standard python list)
# Input:
#   - var: The variable to check for arrayness
# Output:
#   - isArray: True when the function argument is an array
def IsArray(var):
    return isinstance(var, list) or isinstance(var, np.ndarray)

# StringPad is a utility function for fancy string formatting to increase
# readability of output data.
# Input:
#   - name: The prepended piece of text in front of the value argument
#   - value: The value to display
#   - decimals: The number of decimals to display
#   - length: The total length of the resulting stringified value (including
#       decimals). If the displayed number of decimals is lower than the
#       specified number of decimals and length is premitting then the decimals
#       will have '0's appended.
# Output:
#   - result: A formatted string
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
    else:
        stringCenter = ver + ' ' + name

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
              cmap='gnuplot2', contours=None, forceNormMin=None, forceNormMax=None,
              fontsize=15):
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
    axImage.set_xlabel(xLabel, fontsize=fontsize)
    axImage.set_ylabel(yLabel, fontsize=fontsize)
    axImage.grid(True)

    if contours != None:
        ct = axImage.contour(xAxis, yAxis, data, contours, colors='k')
        axImage.clabel(ct)

    # Add the colorbar
    cbb = mpl.colorbar.ColorbarBase(axColorbar, cmap=cmap, norm=norm)
    cbb.set_label(dataLabel, rotation=90, fontsize=fontsize)

def IsAscending(val):
    for i in range(0, len(val) - 1):
        if val[i + 1] < val[i]:
            return False

    return True

def IsDescending(val):
    for i in range(0, len(val) - 1):
        if val[i + 1] > val[i]:
            return False

    return True

def Find1DBisectionAscending(axis, target):
    left = int(0)
    right = int(len(axis) - 1)

    for i in range(0, 128): # An arbitrary limit to guard against improper use
        # Calculate center value and compare
        center = int((left + right) / 2)
        value = axis[center]

        if target < value:
            # Need to look to the left of the center
            if center - left <= 1:
                return left

            right = center
        else:
            # Need to look to the right of the center
            if right - center <= 1:
                return center

            left = center

    raise ValueError("Could not find target, probable improper use of function")

def Find1DBisectionDescending(axis, target):
    left = int(0)
    right = int(len(axis) - 1)

    for i in range(0, 128): # Arbitrary limit
        # Calculate center value and compare
        center = int((left + right) / 2)
        value = axis[center]

        if target > value:
            # Need to look to the left of the center
            if center - left <= 1:
                return left

            right = center
        else:
            # Need to look to the right of the center
            if right - center <= 1:
                return center

            left = center

    raise ValueError("Could not find target, probable improper use of function")

def AdjustSeverity(atmosphere, severity):
    if severity > 0.0:
        return atmosphere[1] + (atmosphere[2] - atmosphere[1]) * severity

    return atmosphere[1] + (atmosphere[1] - atmosphere[0]) * severity

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
        curBias = biasMaps[iMap].getBaseBias()
        print('adjusting', names[iMap], 'bias from', round(curBias, 3),
            'to', round(curBias - amount, 3))

        newBiases[iMap] = curBias - amount
        curBias.reset(newBiases[iMap])

    return newBiases

def Lerp(x1, y1, x2, y2, x):
    return y1 + (y2 - y1) / (x2 - x1) * (x - x1)

# Bilerp will perform bilinear interpolation on the four provided points and
# the two provided coordinates. The z-values are the values that are being
# interpolated. The points should be specified as:
#
#  2 ----- 3    Hence: Point 1 should be at the bottom-left (minimum x and
#  |       |    minimum y-coordinate) and then proceed clockwise. This ordering
#  |       |    will not be checked as the function is called often and
#  1 ----- 4    checking for validity will detract from performance. In case it
#               is not clear: point 1 and 2 are on the left, pointer 3 and 4
# are on the right. Point 1 and 4 are at the bottom and point 2 and 3 are at
# the top.
def Bilerp(xl, xr, yb, yt, z1, z2, z3, z4, x, y):
    # Perform linear interpolation along the x-direction
    zBottom = Lerp(xl, z1, xr, z4, x)
    zTop = Lerp(xl, z2, xr, z3, x)
    return Lerp(yb, zBottom, yt, zTop, y)

def FindXBilerp(xl, xr, yb, yt, z1, z2, z3, z4, zTarget, y):
    zLeft = Lerp(yb, z1, yt, z2, y)
    zRight = Lerp(yb, z4, yt, z3, y)
    return Lerp(zRight, xr, zLeft, xl, zTarget)

def CumulativeSimps(y, x):
    if len(y) != len(x):
        raise ValueError("Expected 'y' and 'x' to be of the same length")

    result = np.zeros(len(y))
    numWhole = int((len(y) - 1)/ 2)

    for i in range(0, numWhole):
        result[2 * i + 1] = result[2 * i] + \
            0.5 * (y[2 * i] + y[2 * i + 1]) * (x[2 * i + 1] - x[2 * i])
        result[2 * i + 2] = result[2 * i] + \
            (x[2 * i + 2] - x[2 * i]) / 6.0 * (
                y[2 * i] + 4.0 * y[2 * i + 1] + y[2 * i + 2])

    if len(y) % 2 == 0:
        # Even number, end with a single trapz section
        result[-1] = result[-2] + 0.5 * (y[-2] + y[-1]) * (x[-1] - x[-2])

    return result

def FormatTime(numSeconds, format):
    # Extract hours, minutes and seconds
    remaining = int(numSeconds)
    seconds = int(numSeconds % 60)
    remaining -= seconds
    minutes = int((remaining / 60) % 60)
    remaining -= minutes * 60
    hours = int(remaining / 3600)

    format = format.lower()

    if format.find('hh') != -1:
        strHours = str(hours)

        while len(strHours) < 2:
            strHours = '0' + strHours

        format = format.replace('hh', strHours)

    if format.find('mm') != -1:
        strMinutes = str(minutes)

        while len(strMinutes) < 2:
            strMinutes = '0' + strMinutes

        format = format.replace('mm', strMinutes)

    if format.find('ss') != -1:
        strSeconds = str(seconds)

        while len(strSeconds) < 2:
            strSeconds = '0' + strSeconds

        format = format.replace('ss', strSeconds)

    return format

def __TestCumulativeSimps__():
    x = np.linspace(0, 4 * np.pi, 1000)
    yInput = - np.sin(x)
    yCorrect = np.cos(x)
    yIntegrate = CumulativeSimps(yInput, x)

    fig = plt.figure()
    ax = fig.add_subplot(121)
    ax.plot(x, yCorrect, 'g')
    ax.plot(x, yIntegrate + 1, 'r')
    ax.grid(True)

    ax = fig.add_subplot(122)
    ax.plot(x, yCorrect - yIntegrate - 1, 'r--')
    ax.grid(True)

#__TestLoadAerodynamicData__()
#__TestCumulativeSimps__()
