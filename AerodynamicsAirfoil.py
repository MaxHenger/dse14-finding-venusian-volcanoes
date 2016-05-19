# -*- coding: utf-8 -*-
"""
Created on Mon May 16 22:52:15 2016

@author: MaxHenger
"""

class AirfoilGenerator:
    def GetUpperSurface(self):
        raise RuntimeError("GetUpperSurface() is called in base class 'AirfoilGenerator'")
        
    def GetLowerSurface(self):
        raise RuntimeError("GetLowerSurface() is called in base class 'AirfoilGenerator'")
        
    def GetCamberSurface(self):
        raise RuntimeError("GetCamberSurface() is called in base calss 'AirfoilGenerator'")
        
    def GetIdentifier(self):
        raise RuntimeError("GetIdentifier() is called in base class 'AirfoilGenerator'")