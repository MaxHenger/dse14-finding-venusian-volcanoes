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