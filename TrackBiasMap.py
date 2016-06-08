# -*- coding: utf-8 -*-
"""
Created on Mon May 30 11:45:23 2016

@author: MaxHenger
"""

import TrackLookup

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as scp_spc

class BiasMap:
    def __init__(self, name, lower, upper, num, biasBase):
        self._name_ = name
        self._axis_ = np.linspace(lower, upper, num)
        self._search_ = None
        self._map_ = np.repeat(biasBase, num)
        self._base_ = biasBase

        if upper > lower:
            self._search_ = TrackLookup.__find1DBisectionAscending__
        else:
            self._search_ = TrackLookup.__find1DBisectionDescending__

    def __call__(self, value):
        index = self._search_(self._axis_, value)
        return TrackLookup.__1DLerp__(self._axis_[index], self._map_[index],
            self._axis_[index + 1], self._map_[index + 1], value)

    def getName(self):
        return self._name_

    def getLower(self):
        return self._axis_[0]

    def getUpper(self):
        return self._axis_[0]

    def getBaseBias(self):
        return self._base_

    def getAxis(self):
        return self._axis_

    def getMap(self):
        return self._map_

    def _erf_(self, biasOffset, location, width):
        return biasOffset * (1 - np.power(scp_spc.erf((self._axis_ - location) / (0.5 * width)), 2.0))

    def reset(self, biasBase):
        self._base_ = biasBase
        self._map_ = np.repeat(biasBase, len(self._axis_))

    def modifyMinimum(self, biasOffset, location, width):
        erfValues = self._erf_(biasOffset, location, width)
        allPositive = True

        for i in range(0, len(self._axis_)):
            if self._axis_[i] <= location:
                self._map_[i] -= biasOffset
            else:
                self._map_[i] -= erfValues[i]

            if self._map_[i] < 0.0:
                allPositive = False

        return allPositive

    def modifyMaximum(self, biasOffset, location, width):
        erfValues = self._erf_(biasOffset, location, width)
        allPositive = True

        for i in range(0, len(self._axis_)):
            if self._axis_[i] >= location:
                self._map_[i] -= biasOffset
            else:
                self._map_[i] -= erfValues[i]

            if self._map_[i] < 0.0:
                allPositive = False

        return allPositive

    def modifyCentered(self, biasOffset, location, width):
        erfValues = self._erf_(biasOffset, location, width)
        allPositive = True

        for i in range(0, len(self._axis_)):
            self._map_[i] -= erfValues[i]

            if self._map_[i] < 0.0:
                allPositive = False

        return allPositive

    def modifyGlobal(self, biasOffset):
        allPositive = True

        for i in range(0, len(self._axis_)):
            self._map_[i] -= biasOffset

            if self._map_[i] < 0.0:
                allPositive = False

        return allPositive

def __TestBiasMap__():
    bm = BiasMap("gamma", 1, 5, 100, 0.9)
    fig = plt.figure()
    ax = fig.add_subplot(131)
    ax.plot(bm.getAxis(), bm.getMap(), 'r')
    print('min =   ', bm.modifyMinimum(0.5, 1.5, 1.1))
    ax.plot(bm.getAxis(), bm.getMap(), 'g')
    print('center =', bm.modifyCentered(0.3, 3.0, 1.0))
    ax.plot(bm.getAxis(), bm.getMap(), 'b')
    print('global =', bm.modifyGlobal(0.5))
    ax.plot(bm.getAxis(), bm.getMap(), 'k')
    ax.grid(True)

    ax = fig.add_subplot(132)
    bm = BiasMap("alpha", -5, 5, 200, 0.95)
    ax.plot(bm.getAxis(), bm.getMap(), 'r')
    print('max =   ', bm.modifyMaximum(0.95, 3.5, 1.0))
    ax.plot(bm.getAxis(), bm.getMap(), 'g')
    print('center =', bm.modifyCentered(0.4, 0, 1.5))
    ax.plot(bm.getAxis(), bm.getMap(), 'b')
    print('center =', bm.modifyCentered(0.5, -2.5, 2.5))
    ax.plot(bm.getAxis(), bm.getMap(), 'k')

    bm = BiasMap("beta", -5, 5, 113, 0.95)
    bm.modifyCentered(0.5, 0, 3.5)
    x = np.linspace(-5, 5, 13)
    y = np.zeros(x.shape)

    for i in range(0, len(x)):
        y[i] = bm(x[i])

    ax = fig.add_subplot(133)
    ax.plot(bm.getAxis(), bm.getMap(), 'g')
    ax.plot(x, y, 'r')

#__TestBiasMap__()
