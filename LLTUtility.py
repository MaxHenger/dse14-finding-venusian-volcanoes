# -*- coding: utf-8 -*-
"""
Created on Wed May 18 22:54:20 2016

@author: MaxHenger
"""

import numpy as np
import matplotlib.pyplot as plt

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

            iLeft = iCenter
        else:
            if iCenter - iLeft <= 1:
                if abs(toFind - array[iLeft]) < abs(toFind - current):
                    return iLeft

                return iCenter

            iRight = iCenter

def stripSpaces(string):
    result = ''

    for i in range(0, len(string)):
        if string[i] != ' ':
            result += string[i]

    return result

def lerp(x1, y1, x2, y2, x):
    return y1 + (y2 - y1) / (x2 - x1) * (x - x1)

def interpolate(x, y, targetX):
    iLeft = 0
    iRight = len(x) - 1

    while True:
        iCenter = int((iLeft + iRight) / 2)
        current = x[iCenter]

        if current <= targetX:
            if iRight - iCenter <= 1:
                return lerp(x[iCenter], y[iCenter], x[iRight], y[iRight], targetX)

            iLeft = iCenter
        else:
            if iCenter - iLeft <= 1:
                return lerp(x[iLeft], y[iLeft], x[iCenter], y[iCenter], targetX)

            iRight = iCenter

def differentiate(x, y):
    dydx = np.zeros(len(x))

    dydx[0] = (y[1] - y[0]) / (x[1] - x[0])
    for i in range(1, len(x) - 1):
        dydx[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
    dydx[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])

    return dydx

def __testLerp__():
    toTest = [
        [0, 1, 0, 2, 0.5, 1],
        [0, 2, 0, 2, 1, 1],
        [0, 3, 0, 1, 1, 0.333333],
        [0, 3, 0, 1, 2, 0.666666],
        [0, 3, 0, 1, 3, 1]
    ]
    
    for i in range(0, len(toTest)):
        ans = lerp(toTest[i][0], toTest[i][2], toTest[i][1], 
                   toTest[i][3], toTest[i][4])
        if abs(ans - toTest[i][5]) > 1e-4:
            print('Error interpolating', toTest[i], 'found', ans)
            
def __testInterpolate__():
    x = np.linspace(0, 200, 11)
    y = np.linspace(0, 200, 11)
    toTest = np.linspace(0, 200, 31)
    
    for i in range(0, 30):
        if abs(toTest[i] - interpolate(x, y, toTest[i])) > 1e-5:
            print('Failed to interpolate:', toTest[i], 
                  'at index', i, 
                  'found', interpolate(x, y, toTest[i]))

def __testDifferentiate__():
    x = np.linspace(0, 2 * np.pi, 80)
    y = np.sin(x)
    plt.plot(x, differentiate(x, y), 'r')
    plt.plot(x, np.cos(x), 'g')
    
def __testClosest__():
    x = np.linspace(0, 10, 11)
    xUnder = np.linspace(-0.1, 9.9, 11)
    xOver = np.linspace(0.1, 10.1, 11)
    
    for i in range(0, 11):
        if findClosestAscending(x, xUnder[i]) != i:
            print("Failed to find closest under ascending", i)
        
        if findClosestUnordered(x, xUnder[i]) != i:
            print("Failed to find closest under unordered", i)
            
        if findClosestAscending(x, xOver[i]) != i:
            print("Failed to find closest over ascending", i)
            
        if findClosestUnordered(x, xOver[i]) != i:
            print("Failed to find closest over unordered", i)
            
#__testLerp__()
#__testInterpolate__()
#__testClosest__()
#__testDifferentiate__()