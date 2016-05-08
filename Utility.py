# -*- coding: utf-8 -*-
"""
Created on Wed May  4 16:57:58 2016

@author: MaxHenger
"""

import numpy as np

def isAscending(array):
    if len(array) == 0:
        raise ValueError("Expected at least one value in 'array'")
        
    prev = array[0]
    
    for i in range(1, len(array)):
        if array[i] < prev:
            return False
            
        prev = array[i]
            
    return True
    
def isDescending(array):
    if len(array) == 0:
        raise ValueError("Expected at least one value in 'array'")
        
    prev = array[0]
    
    for i in range(1, len(array)):
        if array[i] > prev:
            return False
        
        prev = array[i]
        
    return True

def isArray(array):
    return isinstance(array, list) or isinstance(array, np.ndarray)