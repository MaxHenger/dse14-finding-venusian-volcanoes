# -*- coding: utf-8 -*-
"""
Created on Sun May 29 19:25:25 2016

@author: MaxHenger
"""

import numpy as np
import TrackLookup

def LoadAerodynamicData(dataCl, dataCd, tol=1e-15):
    # Load data from file
    dataCl = np.genfromtxt(dataCl, delimiter=';')
    dataCd = np.genfromtxt(dataCd, delimiter=';')
    
    # Check for consistency of dimensions
    if dataCl.shape[1] < 2:
        raise ValueError("Expected at least two columns in Cl data file")
        
    if dataCd.shape[1] < 2:
        raise ValueError("Expected at least two columns in Cd data file")
        
    if dataCl.shape[0] != dataCd.shape[0]:
        raise ValueError("Expected same number of rows in Cl as in Cd file")
        
    # Check for consistency of angles of attack
    for iAlpha in range(0, dataCl.shape[0]):
        if abs(dataCl[iAlpha, 0] - dataCd[iAlpha, 0]) > tol:
            raise ValueError('ClAlpha =', dataCl[iAlpha, 0], '!=',
                             'CdAlpha =', dataCd[iAlpha, 0], 'at index', iAlpha)
            
    # Create the lookup tables and return them
    return TrackLookup.Lookup1D(dataCl[:, 0], dataCl[:, 1]), \
        TrackLookup.Lookup1D(dataCd[:, 0], dataCd[:, 1])
        
def __TestLoadAerodynamicData__():
    Cl, Cd = LoadAerodynamicData("./data/aerodynamicPerformance/Cl.csv",
                                 "./data/aerodynamicPerformance/Cd.csv")
    
    ClPoints = Cl.getPoints()
    CdPoints = Cd.getPoints()
    
    for i in range(0, len(ClPoints[0])):
        print('alpha =', round(ClPoints[0][i], 4),
              ', Cl =', round(ClPoints[1][i], 4),
              ', Cd =', round(CdPoints[1][i], 4))
        
#__TestLoadAerodynamicData__()