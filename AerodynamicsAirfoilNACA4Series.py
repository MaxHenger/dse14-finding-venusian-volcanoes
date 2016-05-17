# -*- coding: utf-8 -*-
"""
Created on Mon May 16 22:56:30 2016

@author: MaxHenger
"""

import numpy as np
import matplotlib.pyplot as plt
import AerodynamicsUtil as aeroutil
import AerodynamicsAirfoil as aerofoil
import AerodynamicsGenerator as aerogen

class AirfoilNACA4Series(aerofoil.AirfoilGenerator):
    def __init__(self, code, chord, close=False):
        # Check the input for validity
        if isinstance(code, str):
            if len(code) != 4:
                raise ValueError("Expected variable 'code' to be a string of " +
                    "exacty 4 characters")
                
            for i in range(0, len(code)):
                if code[i] < '0' or code[i] > '9':
                    raise ValueError("Expected variable 'code' to be a string " +
                        "containing only numbers")
                    
            code = int(code)
        elif isinstance(code, int):
            if code >= 10000:
                raise ValueError("Expected variable 'code' to be smaller than " +
                    "10000, as to describe a NACA 4-digit airfoil")
        else:
            raise ValueError("Expected variable 'code' to be a string or an " +
                "integer")
 
        if aeroutil.isArray(chord):
            if not aeroutil.isAscending(chord):
                raise ValueError("Expected variable 'chord' to be ascending")
        elif isinstance(chord, aerogen.Generator):
            chord = chord.Get()
        else:                
            raise ValueError("Expected variable 'chord' to be a list or " +
                "a subclass of AerodynamicsGenerator.Generator")
        
        # Store the identifying name and the chord positions
        self.identifier = "NACA " + aeroutil.padLeading(str(code), 4, '0')
        self.chord = chord
        self.chordOffset = min(chord)
        self.chordLength = max(chord) - self.chordOffset
        self.chordBase = (self.chord - self.chordOffset) / self.chordLength
        
        # Generate the thickness and the mean chord line
        self.thicknessMaximum = (code % 100) / 100
        self.thicknessHalf = 5 * self.thicknessMaximum * self.chordLength * \
            (
                0.2969 * np.sqrt(self.chordBase) -
                0.1260 * self.chordBase -
                0.3516 * np.power(self.chordBase, 2.0) +
                0.2843 * np.power(self.chordBase, 3.0) -
                0.1015 * np.power(self.chordBase, 4.0)
            )
        
        if (close):
            self.thicknessHalf[-1] = 0.0
            
        self.camberMaximumLocation = ((code - self.thicknessMaximum) % 1000) / 1000.0
        self.camberMaximum = (code - code % 1000) / 100000.0
        self.camber = np.zeros(self.chord.shape)
        self.camberDerivative = np.zeros(self.chord.shape)
        
        for i in range(0, len(self.camber)):
            if (self.chordBase[i] < self.camberMaximumLocation * self.chordLength):
                self.camber[i] = self.camberMaximum * self.chordLength * self.chordBase[i] / \
                    np.power(self.camberMaximumLocation, 2.0) * (
                        2.0 * self.camberMaximumLocation -
                        self.chordBase[i]
                    )
                self.camberDerivative[i] = 2.0 * self.camberMaximum / \
                    np.power(self.camberMaximumLocation, 2.0) * (
                        self.camberMaximumLocation - self.chordBase[i]
                    )
            else:
                self.camber[i] = self.camberMaximum * self.chordLength * (1.0 - self.chordBase[i]) / \
                    np.power(1 - self.camberMaximumLocation, 2.0) * (
                        1.0 + self.chordBase[i] - 2.0 * self.camberMaximumLocation                    
                    )
                self.camberDerivative[i] = 2.0 * self.camberMaximum / \
                    np.power(1 - self.camberMaximumLocation, 2.0) * (
                        self.camberMaximumLocation - self.chordBase[i]                   
                    )
        
        # Construct airfoil from the calculate mean camber line and half thickness
        theta = np.arctan(self.camberDerivative)
        sinTheta = np.sin(theta)
        cosTheta = np.cos(theta)
    
        self.chordUpper = self.chord - (self.thicknessHalf * sinTheta)
        self.chordLower = self.chord + self.thicknessHalf * sinTheta
        self.surfaceUpper = self.camber + self.thicknessHalf * cosTheta
        self.surfaceLower = self.camber - self.thicknessHalf * cosTheta
        
        # Offset the chord such that it starts at the first value the user
        # provided
        self.chordUpper += (self.chordOffset - min(self.chordUpper))
        self.chordLower += (self.chordOffset - min(self.chordLower))
        
    def GetUpperSurface(self):
        return np.transpose([self.chordUpper, self.surfaceUpper])
        
    def GetLowerSurface(self):
        return np.transpose([self.chordLower, self.surfaceLower])
        
    def GetIdentifier(self):
        return self.identifier
        
def __testAirfoilNACA4Series__(code, generator):
    foil = AirfoilNACA4Series(code, generator, True)
    upper = foil.GetUpperSurface()
    lower = foil.GetLowerSurface()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(upper[:, 0], upper[:, 1], 'rx')
    ax.plot(lower[:, 0], lower[:, 1], 'gx')
    ax.axis('equal')
    ax.set_title(foil.GetIdentifier())
    
#__testAirfoilNACA4Series__('2412', aerogen.GeneratorSingleChebyshev(0.0, 1.0, 50))
#__testAirfoilNACA4Series__('4812', aerogen.GeneratorDoubleChebyshev(0.0, 1.0, 50))
#__testAirfoilNACA4Series__('4812', aerogen.GeneratorDoubleChebyshev(1.0, 2.0, 50))
#__testAirfoilNACA4Series__('1236', aerogen.GeneratorLinear(0.0, 1.0, 50))