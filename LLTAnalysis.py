# -*- coding: utf-8 -*-
"""
Created on Thu May 19 00:03:15 2016

@author: MaxHenger
"""

import LLTUtility as lltutil
import LLTAirfoilDatabase as lltdb

import numpy as np
import matplotlib.pyplot as plt

class Wing:
    def __init__(self, airfoil, chordRoot, chordTip, span, numSpan):
        if not isinstance(airfoil, lltdb.Airfoil):
            raise ValueError("Expected variable 'airfoil' to be an instance of the " +
                "LLTAirfoilDatabase.Airfoil class")
            
        if numSpan % 2 == 0:
            numSpan += 1
            
        self.__airfoil__ = airfoil
        self.__chordRoot__ = chordRoot
        self.__chordTip__ = chordTip
        self.__span__ = span
        self.__numSpan__ = numSpan
        
    def analyze(self, alphaMin, alphaMax, numAlpha, Re, numIterations, iterFactor=0.15):
        alpha = np.linspace(alphaMin, alphaMax, numAlpha)
        
        chord = np.append(np.linspace(self.__chordTip__, self.__chordRoot__, 
                                      (self.__numSpan__ - 1) / 2, endpoint=False),
                          np.linspace(self.__chordTip__, self.__chordRoot__, 
                                      (self.__numSpan__ + 1) / 2))
                
        span = np.linspace(-self.__span__ / 2.0, self.__span__ / 2.0, self.__numSpan__)
        
        deltaY = self.__span__ / (self.__numSpan__ - 1)
        gamma = np.zeros([len(alpha), self.__numSpan__])
        
        for i in range(0, len(alpha)):
            curAlpha = alpha[i]
            
            # Take an initial guess at the gamma over velocity distribution
            gammaOverVelocity = np.zeros(self.__numSpan__)
            
            for j in range(0, len(span)):
                gammaOverVelocity[j] = 0.5 * chord[(self.__numSpan__ + 1) / 2] * \
                    np.sqrt(1 - (2 * span[j] / self.__span__)**2) * \
                    self.__airfoil__.GetCl(Re, curAlpha)
                    
            dGammady = lltutil.differentiate(span, gammaOverVelocity)
            print(dGammady)
            
            # Start iterating
            for j in range(0, numIterations):
                # Calculate the induced angle of attack
                alphaInduced = np.zeros(self.__numSpan__)
                
                for k in range(0, self.__numSpan__):
                    value = 0.0
                    
                    for l in range(1, self.__numSpan__, 2):
                        if (abs(span[k] - span[l - 1]) < 1e-5):
                            value += 0.5 * (
                                dGammady[l - 1] / (span[k - 1] - span[l - 1]) +
                                dGammady[l - 1] / (span[k + 1] - span[l - 1])
                            )
                        else:
                            value += dGammady[l - 1] / (span[k] - span[l - 1])
                            
                        if (abs(span[k] - span[l]) < 1e-5):
                            value += 2.0 * (
                                dGammady[l] / (span[k - 1] - span[l]) + \
                                dGammady[l] / (span[k + 1] - span[l])
                            )
                        else:
                            value += 4.0 * dGammady[l] / (span[k] - span[l])
                            
                        if (abs(span[k] - span[l + 1]) < 1e-5):
                            value += 0.5 * (
                                dGammady[l + 1] / (span[k - 2] - span[l + 1]) +
                                dGammady[l + 1] / (span[k - 1] - span[l + 1])
                            )
                        else:
                            value += dGammady[l + 1] / (span[k] - span[l + 1])
                            
                    value *= deltaY / (12.0 * np.pi)
                    alphaInduced[k] = value
                
                # Update gamme over velocity estimate and its derivative
                #print('cur alpha', curAlpha, ', induced', alphaInduced)
                alphaEffective = curAlpha - alphaInduced
                gammaOverVelocityNew = np.zeros(self.__numSpan__)
                
                for k in range(0, self.__numSpan__):
                    gammaOverVelocityNew[k] = 0.5 * chord[k] * self.__airfoil__.GetCl(Re, alphaEffective[k])
                
                gammaOverVelocity = gammaOverVelocity + iterFactor * (gammaOverVelocityNew - gammaOverVelocity)
                dGammady = lltutil.differentiate(span, gammaOverVelocity)
            
            gamma[i] = gammaOverVelocity
            
        fig = plt.figure()
        ax = fig.add_subplot(111)
        lines = []
        desc = []
        
        for i in range(0, len(alpha)):
            if i == 0:
                print(alpha)
                print(gamma[i, :])
                
            curLine, = ax.plot(span, gamma[i, :])
            lines.append(curLine)
            desc.append("a = " + str(alpha[i]))
        
        ax.legend(lines, desc)
        
def __testWing__():
    db = lltdb.Database()
    db.loadFile("airfoilNACA1412.txt")
    w = Wing(db.__airfoils__[0], 2.0, 1.0, 4.0, 100)
    w.analyze(0, 10, 6, 1000000, 200)
    
__testWing__()