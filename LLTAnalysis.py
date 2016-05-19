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
                          np.linspace(self.__chordRoot__, self.__chordTip__, 
                                      (self.__numSpan__ + 1) / 2))
                
        span = np.linspace(-self.__span__ / 2.0, self.__span__ / 2.0, self.__numSpan__)
        
        deltaY = self.__span__ / (self.__numSpan__ - 1)
        gamma = np.zeros([len(alpha), self.__numSpan__])
        
        
        # REMOVE THIS
        initialGuess = np.zeros([len(alpha), 2, self.__numSpan__])
        finalAlpha = np.zeros([len(alpha), self.__numSpan__])
        
        for i in range(0, len(alpha)):
            curAlpha = alpha[i]
            
            # Take an initial guess at the gamma over velocity distribution
            gammaOverVelocity = np.zeros(self.__numSpan__)
            
            for j in range(0, len(span)):
                gammaOverVelocity[j] = 0.5 * chord[int((self.__numSpan__ + 1) / 2)] * \
                    (1 - (2 * span[j] / self.__span__)**2) * \
                    self.__airfoil__.GetCl(Re, curAlpha)
                    
            initialGuess[i, 0] = gammaOverVelocity
                    
            dGammady = lltutil.differentiate(span, gammaOverVelocity)
            
            initialGuess[i, 1] = dGammady
            
            # Start iterating
            for j in range(0, numIterations):
                # Calculate the induced angle of attack
                alphaInduced = np.zeros(self.__numSpan__)
                
                # First calculate all summation terms
                for k in range(0, self.__numSpan__):
                    value = 0.0
                    
                    for l in range(0, self.__numSpan__ - 1):
                        if k == l:
                            value += 2.0 * dGammady[l + 1] / (span[k] - span[l + 1])
                        elif k == l + 1:
                            value += 2.0 * dGammady[l] / (span[k] - span[l])
                        else:
                            value += dGammady[l] / (span[k] - span[l]) + \
                                dGammady[l + 1] / (span[k] - span[l + 1])
                                
                    '''for l in range(0, self.__numSpan__ - 1):
                        if k == l:
                            value += 2.0 * gammaOverVelocity[l + 1] / (span[k] - span[l + 1])
                        elif k == l + 1:
                            value += 2.0 * gammaOverVelocity[l] / (span[k] - span[l])
                        else:
                            value += gammaOverVelocity[l] / (span[k] - span[l]) + \
                                gammaOverVelocity[l + 1] / (span[k] - span[l + 1])'''
                                
                    #value *= 1.0 / (8.0 * np.pi)
                    value *= 180.0 * deltaY / (8.0 * np.pi**2)
                    alphaInduced[k] = value
                
                '''for k in range(0, self.__numSpan__):
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
                            
                    value *= deltaY / (12.0 * np.pi)'''
                
                # Update gamme over velocity estimate and its derivative
                #print('cur alpha', curAlpha, ', induced', alphaInduced
                alphaEffective = curAlpha - alphaInduced
                finalAlpha[i] = alphaInduced
                gammaOverVelocityNew = np.zeros(self.__numSpan__)
                
                for k in range(0, self.__numSpan__):
                    gammaOverVelocityNew[k] = 0.5 * chord[k] * self.__airfoil__.GetCl(Re, alphaEffective[k])
                
                gammaOverVelocity = gammaOverVelocity + iterFactor * (gammaOverVelocityNew - gammaOverVelocity)
                dGammady = lltutil.differentiate(span, gammaOverVelocity)
            
            gamma[i] = gammaOverVelocity
            
        fig = plt.figure()
        ax = fig.add_subplot(141)
        lines = []
        desc = []
        
        for i in range(0, len(alpha)):
            curLine, = ax.plot(span, gamma[i, :])
            lines.append(curLine)
            desc.append("a = " + str(alpha[i]))
        
        ax.legend(lines, desc)
        
        ax = fig.add_subplot(142)
        lines = []
        desc = []
        
        for i in range(0, len(alpha)):
            curLine, = ax.plot(span, initialGuess[i, 0, :])
            lines.append(curLine)
            desc.append('a = ' + str(alpha[i]))
            
        ax.legend(lines, desc)
        
        ax = fig.add_subplot(143)
        lines = []
        desc = []
        
        for i in range(0, len(alpha)):
            curLine, = ax.plot(span, initialGuess[i, 1, :])
            lines.append(curLine)
            desc.append('a = ' +  str(alpha[i]))
            
        ax.legend(lines, desc)
        
        ax = fig.add_subplot(144)
        lines = []
        desc = []
        
        for i in range(0, len(alpha)):
            curLine, = ax.plot(span, finalAlpha[i])
            lines.append(curLine)
            desc.append('a = ' + str(alpha[i]))
            
        ax.legend(lines, desc)
        
def __testWing__():
    db = lltdb.Database()
    db.loadFile("airfoilNACA0012.txt")
    w = Wing(db.__airfoils__[0], 1.0, 0.7, 8.0, 30)
    w.analyze(5, 10, 6, 1000000, 200, 0.05)
    
__testWing__()