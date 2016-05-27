# -*- coding: utf-8 -*-
"""
Created on Tue May 17 10:29:36 2016

@author: MaxHenger
"""

import AerodynamicsUtil as aeroutil

import numpy as np

class Generator:
    def Min(self):
        raise RuntimeError("Min() is called in base class 'Generator'")

    def Max(self):
        raise RuntimeError("Max() is called in base class 'Generator'")

    def Get(self):
        raise RuntimeError("Get() is called in base class 'Generator'")

class GeneratorSingleChebyshev(Generator):
    def __init__(self, minimum, maximum, numPoints):
        self.minimum = minimum
        self.maximum = maximum
        self.values = self.minimum + (self.maximum - self.minimum) * \
            np.cos(np.linspace(np.pi / 2.0, 0, numPoints))

    def Min(self):
        return self.minumum

    def Max(self):
        return self.maximum

    def Get(self):
        return self.values

class GeneratorDoubleChebyshev(Generator):
    def __init__(self, minimum, maximum, numPoints):
        self.minimum = minimum
        self.maximum = maximum
        self.values = self.minimum + (self.maximum - self.minimum) * \
            (np.sin(np.linspace(-np.pi / 2.0, np.pi / 2.0, numPoints)) + 1.0) / 2.0

    def Min(self):
        return self.minimum

    def Max(self):
        return self.maximum

    def Get(self):
        return self.values

class GeneratorLinear(Generator):
    def __init__(self, minimum, maximum, numPoints):
        self.minimum = minimum
        self.maximum = maximum
        self.values = np.linspace(minimum, maximum, numPoints)

    def Min(self):
        return self.minimum

    def Max(self):
        return self.maximum

    def Get(self):
        return self.values

class GeneratorPiecewiseLinear(Generator):
    def __init__(self, values, numPoints):
        if isinstance(numPoints, int):
            numPoints = [numPoints]

        if not aeroutil.isArray(values):
            raise ValueError("Expected variable 'values' to be an array")

        if len(values) != len(numPoints) + 1:
            raise ValueError("Expected the length of the variable 'values' to be " +
                "exactly the length of 'numPoints' + 1")

        self.values = np.linspace(values[0], values[1], numPoints[0], endpoint=False)

        if values[0] < values[1]:
            self.minimum = values[0]
            self.maximum = values[1]
        else:
            self.minimum = values[1]
            self.maximum = values[0]

        for i in range(1, len(numPoints)):
            self.values = np.append(self.values, np.linspace(values[i], values[i + 1], numPoints[i]))

            if values[i + 1] < self.minimum:
                self.minimum = values[i + 1]
            elif values[i + 1] > self.maximum:
                self.maximum = values[i + 1]

    def Min(self):
        return self.minimum

    def Max(self):
        return self.maximum

    def Get(self):
        return self.values

class GeneratorEllipticLinear(Generator):
    def __init__(self, initial, final, numPoints):
        if initial > final:
            self.minimum = final
            self.maximum = initial
        else:
            self.minimum = initial
            self.maximum = final

        self.values = np.zeros(numPoints)

        for i in range(numPoints):
            self.values[i] = final - (final - initial) * \
                np.sqrt(1 - np.power(i / (numPoints - 1), 2.0))

    def Min(self):
        return self.minimum

    def Max(self):
        return self.maximum

    def Get(self):
        return self.values
