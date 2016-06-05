# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 14:40:32 2016

@author: MaxHenger
"""

import TrackStorage
import TrackDiveOptimize
import TrackClimbOptimize

import numpy as np

def CombineTracks(files):
	time = []
	height = []
	vInf = []
	vHor = []
	vVer = []

	lastTime = 0.0

	for i in range(0, len(files)):
		# Load the current file
		data = TrackStorage.DataStorage()
		data.load(files[i])

		# Load the data from the file
		curTime = data.getVariable('time').getValues()
		curHeight = data.getVariable('height').getValues()
		curVInf = data.getVariable('vInf').getValues()
		curVHor = data.getVariable('vHor').getValues()
		curVVer = data.getVariable('vVer').getValues()
		deltaTime = data.getVariable('dt').getValues()

		# Extend the total arrays
		time.extend(curTime + lastTime)
		height.extend(curHeight)
		vInf.extend(curVInf)
		vHor.extend(curVHor)
		vVer.extend(curVVer)

		lastTime = curTime[-1] + deltaTime

	# Store the data
	np.savetxt("JuliusAwesomeFile.csv", np.asarray([time, height, vInf, vHor, vVer]).transpose(),
		'%5.5f', delimiter=';', header='time [s]; height [m]; vInf [m/s]; ' +
		'vHor [m/s]; vVer [m/s]')

CombineTracks(['dive_62000to38000_30to-20_0to0.dat',
			   'climb_38000to62000_40000_-10_0.dat'])
