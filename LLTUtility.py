# -*- coding: utf-8 -*-
"""
Created on Wed May 18 22:54:20 2016

@author: MaxHenger
"""

def findClosestUnordered(array, toFind):
    delta = abs(toFind - array[0])
    best = 0

    for i in range(1, len(array)):
        curDelta = abs(toFind - array[i])

        if curDelta < delta:
            best = i
            delta = curDelta

    return best

def findClosestAscending(array, toFind):
    iLeft = 0
    iRight = len(array) - 1

    while True:
        iCenter = int((iLeft + iRight) / 2)
        current = array[iCenter]

        if current <= toFind:
            if iRight - iCenter <= 1:
                if abs(toFind - array[iRight]) < abs(toFind - current):
                    return iRight

                return iCenter

            iRight = iCenter
        else:
            if iCenter - iLeft <= 1:
                if abs(toFind - array[iLeft]) < abs(toFind - current):
                    return iLeft

                return iCenter

            iLeft = iCenter

def stripSpaces(string):
    result = ''

    for i in range(0, len(string)):
        if string[i] != ' ':
            result += string[i]

    return result

def lerp(x1, y1, x2, y2, x):
    return y1 + (y2 - y1) / (x2 - x1) * x

def interpolate(x, y, targetX):
    iLeft = 0
    iRight = len(x) - 1

    while True:
        iCenter = int((iLeft + iRight)) / 2
        current = x[iCenter]

        if current <= targetX:
            if iRight - iCenter <= 1:
                return lerp(x[iCenter], y[iCenter], x[iRight], y[iRight], targetX)

            iRight = iCenter
        else:
            if iCenter - iLeft <= 1:
                return lerp(x[iLeft], y[iLeft], x[iCenter], y[iCenter], targetX)

            iLeft = iCenter

def differentiate(x, y):
    dydx = np.zeros(len(x))

    dydx[0] = (y[1] - y[0]) / (x[1] - x[0])
    for i in range(1, len(x) - 1):
        dydx[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
    dydx[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])

    return dydx
