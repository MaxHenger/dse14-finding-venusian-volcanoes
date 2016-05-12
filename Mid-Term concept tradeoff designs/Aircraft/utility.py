# -*- coding: utf-8 -*-
"""
Created on Mon May  9 20:05:54 2016

@author: MaxHenger
"""

import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np
import os
import pickle
import ast

class settings:
    def __init__(self):
        # Settable parameters for this specific mission
        self.ratioMass0 = 5.0
        self.ratioMass1 = 0.0
        self.ratioPower0 = 3.0
        self.ratioPower1 = 3.0
        self.ratioTime = 0.3
        self.deltaLatitude = 22.5 * np.pi / 180.0
        
        self.heightLower = 40000
        self.heightUpper = 60000
        
        self.timeLower = 600.0
        self.deltaVUpper = 70
        self.radiusHeatshield = 2.5
        self.massPayload = 100.0
        self.powerPayload = 50.0
        self.specificCapacityBattery = 0.46e6
        self.designLiftCoefficient = 0.4
        self.fuselageRadius = 0.4
        self.minimumDrag = 0.02
        self.oswaldFactor = 0.95
        self.propellorArea = np.pi * self.fuselageRadius ** 2.0 * 3.0
        
        self.efficiencyPower = 0.95
        self.efficiencyCharging = 0.99
        self.efficiencySolarPanel = 0.9 # on top of default efficiency
        
        # Venus' parameters
        self.venusMu = 3.24859e14
        self.venusRadius = 6051800
        self.venusOmega = 2 * np.pi / (243.025 * 24.0 * 60.0 * 60.0)
        self.venusFlux = 2600
        
class dictionaryIO:
    def __init__(self):
        self.__namesAxes__ = None
        self.__valuesAxes__ = None
        self.__resultsLinear__ = None
        self.__multFactors__ = None
        self.__sep__ = ';SEP;SEP;'
        
    def __checkCompoundLength__(self):
        if not (self.__resultsLinear__ is None) and \
                not (self.__namesAxes__ is None) and \
                not (self.__valuesAxes__ is None):
            compoundLength = 1.0
            
            for iAxis in range(0, len(self.__valuesAxes__)):
                compoundLength *= len(self.__valuesAxes__[iAxis])
                
            if len(self.__resultsLinear__) != compoundLength:
                raise ValueError("Linear results do not match compound length of all axes")
    
    def __generateMultFactors__(self):
        self.__multFactors__ = [1] * len(self.__valuesAxes__)
        
        for i in range(1, len(self.__valuesAxes__)):
            iRev = len(self.__valuesAxes__) - 1 - i
            self.__multFactors__[iRev] = self.__multFactors__[iRev + 1] * len(self.__valuesAxes__[iRev + 1])
        
    def setAxes(self, names, values):
        if len(names) != len(values):
            raise ValueError("Expected length of 'names' and 'values' to be the same")
        
        self.__namesAxes__ = names
        self.__valuesAxes__ = values
        
        self.__checkCompoundLength__()
        self.__generateMultFactors__()
        
    def setResults(self, results):
        self.__resultsLinear__ = results
        
        self.__checkCompoundLength__
        
    def save(self, filename):
        fh = open(filename, 'w')
        
        # Write the number of axes
        fh.write(str(len(self.__namesAxes__)) + '\n')
        
        for i in range(0, len(self.__namesAxes__)):
            # Write all axes
            fh.write(self.__namesAxes__[i] + '\n')
            fh.write(str(len(self.__valuesAxes__[i])) + '\n')
            
            fh.write(str(self.__valuesAxes__[i][0]))
            
            for j in range(1, len(self.__valuesAxes__[i])):
                fh.write(self.__sep__ + str(self.__valuesAxes__[i][j]))
                
            fh.write('\n')
            
        # Write the data
        for i in range(0, len(self.__resultsLinear__)):
            first = True
            for key, value in self.__resultsLinear__[i].items():
                # Inefficient seperator writing             
                if first != True:
                    fh.write(self.__sep__)
                
                first = False
                fh.write(key + '=' + str(pickle.dumps(value, protocol=0)))
            
            fh.write('\n')
                
    def load(self, filename):
        self.__namesAxes__ = []
        self.__valuesAxes__ = []
        self.__resultsLinear__ = []
        
        fh = open(filename, 'r')
        
        # Read the number of columns
        line = fh.readline()
        if len(line) == 0 or line[0] == '\n':
            raise ValueError('Expected to find number of axes')

        # Read the columns
        numAxes = int(line)
        numLinear = 1
        
        for i in range(0, numAxes):
            line = fh.readline()
            if len(line) == 0 or line[0] == '\n':
                raise ValueError("Expected to find axis name for axis " + str(i))
                
            self.__namesAxes__.append(line[:-1])
            
            line = fh.readline()
            if len(line) == 0 or line[0] == '\n':
                raise ValueError("Expected to find number of values for axis " + str(i) + ": " + self.__namesAxes__[-1])
            
            numValues = int(line)
            numLinear *= numValues
            line = str(fh.readline())
            if len(line) == 0 or line[0] == '\n':
                raise ValueError("Expected to find values for axis " + str(i) + ": " + self.__namesAxes__[-1])
                
            seperated = line.split(self.__sep__)
            values = [0] * numValues
            
            if numValues != len(seperated):
                raise ValueError("Expected " + str(numValues) + " values for axis " + str(i) + ": " + 
                    self.__namesAxes__[-1] + " but got " + str(len(seperated)) + " values")
            
            for j in range(0, len(seperated)):
                values[j] = float(seperated[j])
            
            self.__valuesAxes__.append(values)
        
        for i in range(0, numLinear):
            line = str(fh.readline())
            if len(line) == 0 or line[0] == '\n':
                raise ValueError("Unexpected EOF at value " + str(i) + " of " + str(numLinear))
                
            seperated = line.split(self.__sep__)
            dictionary = {}
            
            for j in range(0, len(seperated)):
                subSeperated = seperated[j].split('=', 1)
                dictionary[subSeperated[0]] = pickle.loads(ast.literal_eval(subSeperated[1]), encoding='ASCII')
            
            self.__resultsLinear__.append(dictionary)
            
        self.__generateMultFactors__()
        
    def getValue(self, indices):
        if len(indices) != len(self.__valuesAxes__):
            raise ValueError("Expected " + str(len(self.__valuesAxes__)) + " indices, " +
                "received " + len(indices) + " indices")
            
        compoundIndex = 0
        
        for i in range(0, len(indices)):
            compoundIndex += self.__multFactors__[i] * indices[i]
            
        return self.__resultsLinear__[compoundIndex]
        
    def getAxisName(self, indexAxis):
        if indexAxis >= len(self.__namesAxes__):
            raise ValueError("Axis index " + str(indexAxis) +
                " exceeds length " + str(len(self.__namesAxes__)))
            
        return self.__namesAxes__[indexAxis]
        
    def getAxisValues(self, indexAxis):
        if indexAxis >= len(self.__valuesAxes__):
            raise ValueError("Axis index " + str(indexAxis) + 
                " exceeds length " + str(len(self.__valuesAxes__)))
            
        return self.__valuesAxes__[indexAxis]
            
def drawVenus(ax, lat, lon, file):
    # Load image and draw it
    image = img.imread(file)
    ax.imshow(image, extent=[-180.0, 180.0, -90, 90])
    
    # Overlay specified latitude and longitude grid
    ax.set_xticks(lon)
    ax.set_yticks(lat)
    ax.grid(True, color=[1.0, 1.0, 1.0])
    
def __testDictionaryIO__():
    test = dictionaryIO()
    test.setAxes(['a', 'b'], [[0, 1, 2], [25, 50, 75]])
    test.setResults([{'a': 2, 'b': 3}, {'a': 2, 'b': 4}, {'a': 2, 'b': 5},
                     {'a': 3, 'b': 3}, {'a': 3, 'b': 4}, {'a': 3, 'b': 5},
                     {'a': 4, 'b': 3}, {'a': 4, 'b': 4}, {'a': 4, 'b': 5}])
    
    test.save('testFile.txt')
    print(' --- before reloading:')
    for i in range(0, 2):
        print('axis 1 name   =', test.getAxisName(0))
        print('axis 1 values =', test.getAxisValues(0))
        print('axis 2 name   =', test.getAxisName(1))
        print('axis 2 values =', test.getAxisValues(1))
        
        for i in range(0, len(test.getAxisValues(0))):
            for j in range(0, len(test.getAxisValues(1))):
                val = test.getValue([i, j])
                print('[', i, ',', j, '] =', val)
           
        if i == 0:
            test = dictionaryIO()
            test.load('testFile.txt')
            print(' --- after reloading:')
            
#__testDictionaryIO__()