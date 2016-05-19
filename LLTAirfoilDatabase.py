# -*- coding: utf-8 -*-
"""
Created on Wed May 18 22:52:10 2016

@author: MaxHenger
"""

import LLTUtility as lltutil

import matplotlib.pyplot as plt
import numpy as np

class Airfoil:
    def __init__(self):
        self.__name__ = ''
        self.__Re__ = []
        self.__alpha__ = []
        self.__Cl__ = []
        self.__Cd__ = []
        self.__Cm__ = []
        
    def SetName(self, name):
        self.__name__ = name
        
    def AddSeries(self, Re, alpha, Cl, Cd, Cm):
        self.__Re__.append(Re)
        self.__alpha__.append(alpha)
        self.__Cl__.append(Cl)
        self.__Cd__.append(Cd)
        self.__Cm__.append(Cm)
    
    def GetName(self):
        return self.__name__
        
    def GetAlphaSequence(self, Re):
        return self.__alpha__[lltutil.findClosestUnordered(self.__Re__, Re)]
        
    def GetClSequence(self, Re):
        return self.__Cl__[lltutil.findClosestUnordered(self.__Re__, Re)]
        
    def GetCl(self, Re, alpha):
        index = lltutil.findClosestUnordered(self.__Re__, Re)
        return lltutil.interpolate(self.__alpha__[index], self.__Cl__[index], alpha)
        
    def GetCdSequence(self, Re):
        return self.__Cd__[lltutil.findClosestUnordered(self.__Re__, Re)]
        
    def GetCd(self, Re, alpha):
        index = lltutil.findClosestUnordered(self.__Re__, Re)
        return lltutil.interpolate(self.__alpha__[index], self.__Cd__[index], alpha)
        
    def GetCmSequence(self, Re):
        return self.__Cm__[lltutil.findClosestUnordered(self.__Re__, Re)]
        
    def GetCm(self, Re, alpha):
        index = lltutil.findClosestUnordered(self.__Re__, Re)
        return lltutil.interpolate(self.__alpha__[index], self.__Cm__[index], alpha)
        
class Database:
    def __init__(self):
        self.__airfoils__ = []
        
    def eof(self, line):
        return len(line) == 0
        
    def loadFile(self, filename):
        # Start reading file
        fh = open(filename, 'r')
        
        airfoil = Airfoil()
        first = True
        loadingValues = False
        Re = 0
        alpha = []
        Cl = []
        Cd = []
        Cm = []
        
        while True:
            # Retrieve line from the file and check contents
            line = fh.readline()
            
            if self.eof(line):
                if first != True:
                    airfoil.AddSeries(Re, alpha, Cl, Cd, Cm)
                    
                break
            
            if len(line) == 1 and line == '\n':
                continue
            
            if line.find("XFOIL") == 0:
                # indicates start of new airfoil, store the current data and
                # clean out the variables
                loadingValues = False
                if first != True:
                    airfoil.AddSeries(Re, alpha, Cl, Cd, Cm)
                    Re = 0
                    alpha = []
                    Cl = []
                    Cd = []
                    Cm = []
                    
                first = False
            elif loadingValues == True:
                loadedValues = list(filter(None, line.split(" ")))
                
                if len(loadedValues) != 7:
                    raise ValueError("Expected 7 columns when loading values on line", line)
                
                alpha.append(float(loadedValues[0]))
                Cl.append(float(loadedValues[1]))
                Cd.append(float(loadedValues[2]))
                Cm.append(float(loadedValues[4]))
            elif line.find("Calculated polar for:") == 0:
                airfoil.SetName(line[len("Calculated polar for: "):-1])
            elif line.find("Mach") == 0:
                iRe = line.find("Re =")
                iNCrit = line.find("Ncrit =")
                
                if iRe == -1 or iNCrit == -1:
                    raise ValueError("Expected 'mach' and 'Ncrit' on a line containing 'mach'")
                
                Re = float(lltutil.stripSpaces(line[iRe + len("Re ="):iNCrit]))
            elif line.find("---") == 0:
                loadingValues = True
        
        self.__airfoils__.append(airfoil)
        
def __testAirfoil__():
    db = Database()
    db.loadFile("airfoilNACA1412.txt")
    foil = db.__airfoils__[0]
    
    alpha = np.linspace(-20, 20, 100)
    Cl = np.zeros(alpha.shape)
    
    for i in range(0, len(alpha)):
        Cl[i] = foil.GetCl(1000000, alpha[i])
    
    plt.plot(alpha, Cl)

#__testAirfoil__()