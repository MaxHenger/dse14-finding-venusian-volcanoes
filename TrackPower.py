# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 15:43:57 2016

@author: MaxHenger
"""

import numpy as np

# Old power estimating code
#class TrackPower:
#	def __init__(self, settings, atmosphere):
#		self.settings = settings
#		self.atmosphere = atmosphere
#		self.constantIndirect = 0.7
#		self.constantUniform = 0.7
#		self.boundUpper = 67500
#		self.boundCenter = 55000
#		self.boundBottom = 30000
#		self.boundDelta = self.boundUpper - self.boundCenter
#
#	# getPowerEfficiency
#	def getPowerEfficiency(self, height, latitude, longitude, alpha, gamma):
#		efficiency = self.atmosphere.solarEfficiencyBoundless(height, latitude,
#			longitude, includeZenith=False)
#
#		# Modify the efficiency with a direct part and an indirect part based on
#		# a bullshit interpolation scheme
#		if height > self.boundUpper:
#			# Do not apply any modifications: Direct and indirect sunlight
#			base = efficiency * np.cos(latitude) * np.cos(longitude + alpha + gamma)
#			return base, base * self.constantIndirect, 0
#		elif height > self.boundCenter:
#			# Direct and uniform lighting
#			baseDirect = efficiency * np.cos(latitude) * np.cos(longitude + alpha +
#				gamma) * (height - self.boundCenter) / self.boundDelta
#			return baseDirect, baseDirect * self.constantIndirect, \
#				efficiency * self.constantUniform * (self.boundUpper - height) / \
#				self.boundDelta
#		elif height > self.boundBottom:
#			# Only uniform lighting
#			return 0, 0, efficiency * self.constantUniform
#
#		# Only direct sunlight on the top part
#		return efficiency * np.cos(latitude) * np.cos(longitude + alpha +
#			gamma), 0, 0

class TrackPower:
    def __init__(self, settings, atmosphere):
        self.settings = settings
        self.atmosphere = atmosphere
        
    def getPowerEfficiency(self, height, latitude, longitude, alpha, gamma):
        # Get the solar efficiency as predicted by the NASA paper and 
        # compensate for the predicted efficiency of the CPV Point Focus Solar
        # Cells.
        efficiency = np.asarray(self.atmosphere.solarEfficiencyBoundless(height, latitude,
            longitude, includeZenith=False)) / 29 * 35
        
        # The 0.77 figure is taken from the "Measurements of Sunlight Flux on 
        # Venus paper by Tomasko et. al.
        factor = np.cos(latitude) * np.cos(longitude + alpha + gamma)
        efficiency *= factor
        return efficiency, efficiency * 0.77, 0