# -*- coding: utf-8 -*-
"""
Created on Mon May 16 23:37:47 2016

@author: MaxHenger
"""

import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def padLeading(msg, length, character):
    return (character * (length - len(msg))) + msg

def isArray(variable):
    return isinstance(variable, list) or isinstance(variable, np.ndarray)

def isAscending(data):
    for i in range(0, len(data) - 1):
        if (data[i] > data[i + 1]):
            return False

    return True;

def set3DAxesEqual(ax, x, y, z):
    rangeX = max(x) - min(x)
    rangeY = max(y) - min(y)
    rangeZ = max(z) - min(z)
    maxRange = max(rangeX, rangeY, rangeZ)
    centerX = np.average(x)
    centerY = np.average(y)
    centerZ = np.average(z)

    ax.set_xlim(centerX - maxRange * 0.5, centerX + maxRange * 0.5)
    ax.set_ylim(centerY - maxRange * 0.5, centerY + maxRange * 0.5)
    ax.set_zlim(centerZ - maxRange * 0.5, centerZ + maxRange * 0.5)

def plot3DColors(ax, x, y, z, indices, colors):
    collection = np.zeros([len(indices), 3, 3])

    for iIndex in range(0, len(indices)):
        for iSubIndex in range(0, 3):
            collection[iIndex, iSubIndex, 0] = x[int(indices[iIndex][iSubIndex])]
            collection[iIndex, iSubIndex, 1] = y[int(indices[iIndex][iSubIndex])]
            collection[iIndex, iSubIndex, 2] = z[int(indices[iIndex][iSubIndex])]

    toPlot = Poly3DCollection(collection, facecolors=colors, edgecolors='k')
    ax.add_collection(toPlot)
