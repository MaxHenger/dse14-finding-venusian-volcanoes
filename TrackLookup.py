# -*- coding: utf-8 -*-
"""
Created on Fri May 27 12:01:24 2016

TrackLookup contains several generalised definitions of lookup tables. The
classes defined in this file assume that you know what you're doing. They will
not guard against improper use everywhere, although they will at some points in
the code.

@author: MaxHenger
"""

import TrackCommon

def __1DNearestNeighbour__(axis1, result1, axis2, result2, target):
    if (target - axis1 < axis2 - target):
        return result1

    return result2

def __1DLerp__(axis1, result1, axis2, result2, target):
    return result1 + (result2 - result1) / (axis2 - axis1) * (target - axis1)

def __derivative1__(axis, value):
    # Perform some checks on the input
    if len(axis) != len(value):
        raise ValueError("Expected the length of 'axis' to be equal to 'value'")

    if len(axis) < 2:
        raise ValueError("Expected at least two points")

    # Preallocate the array and calculate the first and final value (with the
    # 1-order forward and backward equations)
    result = np.zeros([len(axis)])
    result[0] = (value[1] - value[0]) / (axis[1] - axis[0])
    result[-1] = (value[-1] - value[-2]) / (axis[-1] - axis[-2])

    for i in range(1, len(axis) - 1):
        result[i] = (value[i + 1] - value[i - 1]) / (axis[i + 1] - axis[i - 1])

    return result

def __derivative2__(axis, value):
    # Perform some checks on the input
    if len(axis) != len(value):
        raise ValueError("Expected the length of 'axis' to be equal to 'value'")

    if len(axis) < 3:
        raise ValueError("Expected at least three points")

    # Preallocate array and calculate the first and second initial and final
    # values (using the 2-point forward and 1-point central difference)
    result = np.zeros([len(axis)])
    result[0] = (-3 * value[0] + 4 * value[1] - value[2]) / (axis[2] - axis[0])
    result[1] = (value[2] - value[0]) / (axis[2] - axis[0])
    result[-1] = (value[-1] - value[-3]) / (axis[-1] - axis[-3])
    result[-2] = (value[-3] - 4 * value[-2] + 3 * value[-1]) / (axis[-1] - axis[-3])

    for i in range(2, len(axis) - 2):
        result[i] = (value[i - 2] - 8 * value[i - 1] + 8 * value[i + 1] - value[i + 2]) / \
            (3 * (axis[i + 2] - axis[i - 2]))

    return result

class Lookup1D:
    def __init__(self, axis, result, algorithm="bisection", interpolation="linear"):
        self.__axis__ = axis
        self.__result__ = result

        algorithm = algorithm.lower()
        interpolation = interpolation.lower()

        if algorithm == "bisection":
            if TrackCommon.IsAscending(self.__axis__):
                self.__search__ = TrackCommon.Find1DBisectionAscending
            elif TrackCommon.IsDescending(self.__axis__):
                self.__search__ = TrackCommon.Find1DBisectionDescending
            else:
                raise ValueError("Data should be purely ascending or purely descending")

            if interpolation == "linear":
                self.__interpolate__ = __1DLerp__
            elif interpolation == "nearestneighbour" or interpolation == "nearest":
                self.__interpolate__ = __1DNearestNeighbour__
            else:
                raise ValueError("Unkown interpolation method")
        else:
            raise ValueError("Unknown algorithm")

    def __call__(self, target, *args):
        return self.find(target)

    def find(self, target, *args):
        i = self.__search__(self.__axis__, target)
        return self.__interpolate__(self.__axis__[i], self.__result__[i],
                                    self.__axis__[i + 1], self.__result__[i + 1],
                                    target)

    def getPoints(self):
        # Make sure we're not returning a reference to the internal variables
        xAxis = []
        xAxis.extend(self.__axis__)
        yAxis = []
        yAxis.extend(self.__result__)
        return xAxis, yAxis

    def getDerivative(self, algorithm='2central'):
        # A utility function that returns the derivative of the current lookup
        # map. The used algorithm is not very accurate. But should do for now.
        algorithm = algorithm.lower()

        if algorithm == '1central':
            return Lookup1D(self.__axis__, __derivative1__(self.__axis__, self.__result__))
        elif algorithm == '2central':
            return Lookup1D(self.__axis__, __derivative2__(self.__axis__, self.__result__))
        raise ValueError("Unknown algorithm")

class LookupSegmented1D:
    class Basket:
        def __init__(self, axis, result, vMin, vMax, search, interpolate):
            self.__axis__ = axis
            self.__result__ = result
            self.__vMin__ = vMin
            self.__vMax__ = vMax
            self.__search__ = search
            self.__interpolate__ = interpolate

        def contains(self, target):
            return self.__vMin__ <= target and self.__vMax__ >= target

        def find(self, target):
            i = self.__search__(self.__axis__, target)
            return self.__interpolate__(self.__axis__[i], self.__result__[i],
                                        self.__axis__[i + 1], self.__result__[i + 1],
                                        target)

    def __init__(self, axis, result, algorithm="bisection", interpolation="linear"):
        # Preallocate data and clean input
        self.__baskets__ = []
        algorithm = algorithm.lower()
        interpolation = interpolation.lower()

        # Determine search and interpolation algorithm to use
        search = None
        interpolate = None

        if algorithm == "bisection":
            if interpolation == "linear":
                interpolate = __1DLerp__
            elif interpolation == "nearestneighbour" or interpolation == "nearest":
                interpolate = __1DNearestNeighbour__
            else:
                raise ValueError("Unknown interpolation method")
        else:
            raise ValueError("Unknown algorithm")

        # Start generating baskets
        ascending = (axis[1] - axis[0]) > 0
        startIndex = 0

        for i in range(2, len(axis)):
            curAscending = (axis[i] - axis[i - 1]) > 0

            if curAscending != ascending:
                # Construct a new basket. First figure out what the correct
                # searching and interpolation algorithms are
                if ascending:
                    search = TrackCommon.Find1DBisectionAscending
                else:
                    search = TrackCommon.Find1DBisectionDescending

                curAxis = axis[startIndex:i]

                self.__baskets__.append(self.Basket(curAxis, result[startIndex:i],
                                                    min(curAxis), max(curAxis),
                                                    search, interpolate))

                # Update variables for storing the next basket
                ascending = curAscending
                startIndex = i - 1

        # Store the final basket
        if ascending:
            search = TrackCommon.Find1DBisectionAscending
        else:
            search = TrackCommon.Find1DBisectionDescending

        curAxis = axis[startIndex:]
        self.__baskets__.append(self.Basket(curAxis, result[startIndex:],
                                            min(curAxis), max(curAxis),
                                            search, interpolate))

    def __call__(self, target, *args):
        return self.find(target, args)
        
    def find(self, target, *args):
        # Find all instances where the provided target is in a basket
        results = []

        for i in range(0, len(self.__baskets__)):
            if self.__baskets__[i].contains(target):
                results.append(self.__baskets__[i].find(target))

        return results

    def getPoints(self, target):
        if len(self.__baskets__) == 0:
            # No baskets, return nothing
            return [], []

        # Set the points coming from the first basket
        xAxis = []
        yAxis = []

        xAxis.extend(self.__baskets__[0].__axis__)
        yAxis.extend(self.__baskets__[0].__result__)

        # Set the remaining points, never add the first one: it is a copy from
        # the previous basket
        for i in range(1, len(self.__baskets__)):
            xAxis.extend(self.__baskets__[i].__axis__)
            yAxis.extend(self.__baskets__[i].__result__)

        return xAxis, yAxis

class Lookup2D:
    def __init__(self, axisX, axisY, result, algorithm="bisection", interpolation="bilinear"):
        self.__axisX__ = np.asarray(axisX)
        self.__axisY__ = np.asarray(axisY)
        self.__result__ = np.asarray(result)
        
        algorithm = algorithm.lower()
        interpolation = interpolation.lower()
        
        if algorithm == "bisection":
            # Determine wether the x-axis and the y-axis are ascending or
            # descending, throw an error if none of either
            if TrackCommon.IsAscending(self.__axisX__):
                self.__searchX__ = TrackCommon.Find1DBisectionAscending
            elif TrackCommon.IsDescending(self.__axisX__):
                self.__searchX__ = TrackCommon.Find1DBisectionDescending
            else:
                raise ValueError("x-axis should be purely ascending or descending")
                
            if TrackCommon.IsAscending(self.__axisY__):
                self.__searchY__ = TrackCommon.Find1DBisectionAscending
            elif TrackCommon.IsDescending(self.__axisY__):
                self.__searchY__ = TrackCommon.Find1DBisectionDescending
            else:
                raise ValueError("y-axis should be purely ascending or descending")
                
            if interpolation == "bilinear":
                self.__interpolate__ = TrackCommon.Bilerp
            else:
                raise ValueError("Unknown interpolation method")
        else:
            raise ValueError("Unknown algorithm")
            
    def __call__(self, x, y):
        return self.find(x, y)
        
    def find(self, x, y):
        iX = self.__searchX__(self.__axisX__, x)
        iY = self.__searchY__(self.__axisY__, y)
        
        if iX == len(self.__axisX__) - 1:
            # x-coordinate is on the upper edge
            if iY == len(self.__axisY__) - 1:
                # y-coordinate is on the upper edge
                return self.__result__[iX, iY]

            return self.__interpolate__(self.__axisX__[iX - 1], self.__axisX__[iX],
                self.__axisY__[iY], self.__axisY__[iY + 1], self.__result__[iX - 1, iY],
                self.__result__[iX - 1, iY + 1], self.__result__[iX, iY + 1],
                self.__result__[iX, iY], x, y)
        
        if iY == len(self.__axisY__) - 1:
            # Only the y-coordinate is on the upper edge
            return self.__interpolate__(self.__axisX__[iX], self.__axisX__[iX + 1],
                self.__axisY__[iY - 1], self.__axisY__[iY], self.__result__[iX, iY - 1],
                self.__result__[iX, iY], self.__result__[iX + 1, iY], 
                self.__result__[iX + 1, iY - 1], x, y)
        
        # Interpolate somewhere away from the edges
        return self.__interpolate__(self.__axisX__[iX], self.__axisX__[iX + 1],
            self.__axisY__[iY], self.__axisY__[iY + 1], self.__result__[iX, iY],
            self.__result__[iX, iY + 1], self.__result__[iX + 1, iY + 1],
            self.__result__[iX + 1, iY], x, y)
        
    def getPoints(self):
        # Do not return a reference to the internal variables
        xAxis = []
        xAxis.extend(self.__axisX__)
        yAxis = []
        yAxis.extend(self.__axisY__)
        result = []
        result.extend(self.__result__)
        
        return np.asarray(xAxis), np.asarray(yAxis), np.asarray(result)
        
    def getDerivative(self, algorithm='2central', axis='x'):
        algorithm = algorithm.lower()
        axis = axis.lower()
        result = [[0] * len(self.__axisY__)] * len(self.__axisX__)
        
        # Get derivative-generating function
        derivative = None
        
        if algorithm == '1central':
            derivative = __derivative1__
        elif algorithm == '2central':
            derivative = __derivative2__
        else:
            raise ValueError("Unknown algorithm")
        
        # Generate the derivatives
        if axis == 'x':
            # Go through all y-values
            for iY in range(0, len(self.__axisY__)):
                # Assemble values to take the derivative
                values = [0] * len(self.__axisX__)
                
                for iX in range(0, len(self.__axisX__)):
                    values[iX] = self.__result__[iX, iY]
        
                # Take the derivative and store it
                values = derivative(self.__axisX__, values)
                
                for iX in range(0, len(self.__axisX__)):
                    result[iX][iY] = values[iX]
        elif axis == 'y':
            for iX in range(0, len(self.__axisX__)):
                result[iX] = derivative(self.__axisY__, self.__result__[iX])
        
        return Lookup2D(self.__axisX__, self.__axisY__, result)

class Lookup2DReverse:
    def __init__(self, axisX, axisY, result):
        self.__axisX__ = axisX
        self.__axisY__ = axisY
        
        if TrackCommon.IsAscending(self.__axisX__):
            self.__searchX__ = TrackCommon.Find1DBisectionAscending
        elif TrackCommon.IsDescending(self.__axisX__):
            self.__searchX__ = TrackCommon.Find1DBisectionDescending
        else:
            raise ValueError("x-axis should be purely ascending or descending")
            
        if TrackCommon.IsAscending(self.__axisY__):
            self.__searchY__ = TrackCommon.Find1DBisectionAscending
        elif TrackCommon.IsDescending(self.__axisY__):
            self.__searchY__ = TrackCommon.Find1DBisectionDescending
        else:
            raise ValueError("y-axis should be purely ascending or descending")
            
        self.__result__ = result
        
    def __call__(self, target, y):
        return self.find(target, y)
        
    def find(self, target, y):
        # Make a linear interpolation of the values at the given y-value
        #print('looking for:', round(target, 4), 'at', round(y, 4))
        iY = self.__searchY__(self.__axisY__, y)
        prev = TrackCommon.Lerp(self.__axisY__[iY], self.__result__[0, iY],
            self.__axisY__[iY + 1], self.__result__[0, iY + 1], y)
        
        found = []
        
        for iX in range(1, len(self.__axisX__)):
            new = TrackCommon.Lerp(self.__axisY__[iY], self.__result__[iX, iY],
                self.__axisY__[iY + 1], self.__result__[iX, iY + 1], y)
            
            if (prev <= target and new > target) or (prev > target and new < target):
                found.append(TrackCommon.Lerp(prev, self.__axisX__[iX - 1], new,
                    self.__axisX__[iX], target))
                    
            prev = new
                    
        return found
            
# Do some simple testing
import numpy as np

def __testLookup1D__():
    xPossible = [
        np.linspace(0, 10, 10),
        np.linspace(10, 20, 10),
        np.linspace(10, 0, 10),
        np.linspace(20, 10, 10),
        np.linspace(10, -10, 20)
    ]

    xCheck = [
        np.linspace(0, 10, 100),
        np.linspace(10, 20, 100),
        np.linspace(10, 0, 100),
        np.linspace(20, 10, 100),
        np.linspace(10, -10, 200)
    ]

    if len(xCheck) != len(xPossible):
        raise ValueError("Invalid test definition")

    for j in range(0, len(xPossible)):
        lut = Lookup1D(xPossible[j], xPossible[j])
        err = False

        for i in range(0, len(xCheck[j])):
            found = lut.find(xCheck[j][i])
            if abs(xCheck[j][i] - found) > 1e-8:
                err = True
                print("at i =", i, ", expected =", xCheck[j][i], ", received =", found)

        if err:
            raise ValueError("at j =", j, "the test failed")

#__testLookup1D__()

def __testLookupSegmented1D__():
    xPossible = [
        np.append(np.linspace(0, 10, 10), np.linspace(10, 0, 10))
    ]

    yPossible = [
        np.append(np.linspace(0, 10, 10), np.linspace(10, 20, 10))
    ]

    xCheck = [
        np.linspace(0, 10, 100)
    ]

    if len(xCheck) != len(xPossible):
        raise ValueError("Invalid test definition")

    for j in range(0, len(xPossible)):
        lut = LookupSegmented1D(xPossible[j], yPossible[j])
        err = False

        for i in range(0, len(xCheck[j])):
            found = lut.find(xCheck[j][i])
            if len(found) == 1:
                if abs(xCheck[j][i] - found[0]) > 1e-8:
                    err = True
                    print('at i =', i, ', expected =', xCheck[j][i], ", received =", found[0])
            elif len(found) == 2:
                if abs(xCheck[j][i] - found[0]) > 1e-8:
                    err = True
                    print('at i =', i, '(0), expected =', xCheck[j][i], ", received =", found[0])

                if abs(20.0 - xCheck[j][i] - found[1]) > 1e-8:
                    err = True
                    print('at i =', i, '(1), expected =', 20 - xCheck[j][i], ', received =', found[1])
            else:
                print('at i =', i, len(found), 'values were found, expected 1 or 2')

        if err == True:
            raise ValueError("at j =", j, "the test failed")

#__testLookupSegmented1D__()

def __testLookup1DDerivative__():
    x = np.linspace(0, 4 * np.pi, 200)
    y = np.sin(x)
    yDeriv = np.cos(x)

    lookup = Lookup1D(x, y)
    lookupDeriv = lookup.getDerivative('1central')

    for i in range(0, len(x)):
        if abs(lookupDeriv(x[i]) - yDeriv[i]) > 1e-3:
            raise ValueError("at i =", i, "received", lookupDeriv(x[i]),
                "but expected", yDeriv[i])

    lookup = Lookup1D(x, y)
    lookupDeriv = lookup.getDerivative('2central')

    for i in range(0, len(x)):
        if abs(lookupDeriv(x[i]) - yDeriv[i]) > 1e-2:
            raise ValueError("at i =", i, "received", lookupDeriv(x[i]),
                "but expected", yDeriv[i])

#__testLookup1DDerivative__()

#import matplotlib.pyplot as plt

def __testLookup2D__():
    x = np.linspace(0, 5.0, 10)
    y = np.linspace(0, 5.0, 10)
    z = np.zeros([len(x), len(y)])
    
    for iX in range(0, len(x)):
        for iY in range(0, len(y)):
            z[iX, iY] = x[iX] * y[iY]**1.1
    
    lookup = Lookup2D(x, y, z)
    
    numReinterpolate = 100
    zExpected = np.zeros([numReinterpolate, numReinterpolate])
    zLookup = np.zeros(zExpected.shape)
    xNew = np.linspace(0, 5.0, numReinterpolate)
    yNew = np.linspace(0, 5.0, numReinterpolate)
    
    for iX in range(0, len(xNew)):
        for iY in range(0, len(yNew)):
            zExpected[iX, iY] = xNew[iX] * yNew[iY] ** 1.1
            zLookup[iX, iY] = lookup(xNew[iX], yNew[iY])
            
            if abs(zExpected[iX, iY] - zLookup[iX, iY]) > 1e-10:
                print("At", xNew[iX], ",", yNew[iY], "found", zLookup[iX, iY],
                    "but expected", zExpected[iX, iY])
                
    boundsL, boundsR = TrackCommon.ImageAxes(0.0, 1.0, 0.0, 1.0)
    
    fig = plt.figure()
    axL = fig.add_axes(boundsL)
    axR = fig.add_axes(boundsR)
    TrackCommon.PlotImage(fig, axL, axR, x, 'x', y, 'y', zExpected - zLookup, 'z')
                
#__testLookup2D__()

def __TestLookupReverse2D__():
    x = np.linspace(0, 5.0, 10)
    y = np.linspace(0, 5.0, 10)
    z = np.zeros([len(x), len(y)])
    offset = 0.25
    
    for iX in range(0, len(x)):
        for iY in range(0, len(y)):
            z[iX, iY] = x[iX] * y[iY] * 3.0

    scatterOriginalX = []
    scatterOriginalY = []
    scatterNewX = []
    scatterNewY = []

    lookup = Lookup2DReverse(x, y, z)

    for iX in range(0, len(x) - 1):
        for iY in range(0, len(y) - 1):
            toFind = (z[iX, iY] + z[iX, iY + 1]) / 2.0 + \
                     (z[iX + 1, iY] - z[iX, iY]) * offset
                      
            scatterOriginalX.append(x[iX])
            scatterOriginalY.append(y[iY])
            halfY = (y[iY] + y[iY + 1]) / 2.0
            scatterNewX.append(lookup(toFind, halfY))
            scatterNewY.append(halfY)