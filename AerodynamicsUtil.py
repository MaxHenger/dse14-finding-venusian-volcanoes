# -*- coding: utf-8 -*-
"""
Created on Mon May 16 23:37:47 2016

@author: MaxHenger
"""

import numpy as np

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