# -*- coding: utf-8 -*-
"""
Created on Fri May 27 12:01:24 2016

TrackLookup contains several generalised definitions of lookup tables. The
classes defined in this file assume that you know what you're doing. They will
not guard against improper use everywhere, although they will at some points in
the code.

@author: MaxHenger
"""

def isAscending(val):
    for i in range(0, len(val) - 1):
        if val[i + 1] < val[i]:
            return False

    return True

def isDescending(val):
    for i in range(0, len(val) - 1):
        if val[i + 1] > val[i]:
            return False

    return True

def __find1DBisectionAscending__(axis, target):
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

def __find1DBisectionDescending__(axis, target):
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
    def __init__(self, axis, result, ordered=True, algorithm="bisection", interpolation="linear"):
        self.__axis__ = axis
        self.__result__ = result

        algorithm = algorithm.lower()
        interpolation = interpolation.lower()

        if algorithm == "bisection":
            if isAscending(self.__axis__):
                self.__search__ = __find1DBisectionAscending__
            elif isDescending(self.__axis__):
                self.__search__ = __find1DBisectionDescending__
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

    def __call__(self, target):
        return self.find(target)

    def find(self, target):
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
                    search = __find1DBisectionAscending__
                else:
                    search = __find1DBisectionDescending__

                curAxis = axis[startIndex:i]

                self.__baskets__.append(self.Basket(curAxis, result[startIndex:i],
                                                    min(curAxis), max(curAxis),
                                                    search, interpolate))

                # Update variables for storing the next basket
                ascending = curAscending
                startIndex = i - 1

        # Store the final basket
        if ascending:
            search = __find1DBisectionAscending__
        else:
            search = __find1DBisectionDescending__

        curAxis = axis[startIndex:]
        self.__baskets__.append(self.Basket(curAxis, result[startIndex:],
                                            min(curAxis), max(curAxis),
                                            search, interpolate))

    def __call__(self, target):
        return self.find(target)
        
    def find(self, target):
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

__testLookup1DDerivative__()
