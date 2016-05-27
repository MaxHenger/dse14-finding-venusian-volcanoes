# -*- coding: utf-8 -*-
"""
Created on Thu May 19 13:13:10 2016

@author: MaxHenger
"""

import AerodynamicsWing as aerowing
import AerodynamicsGenerator as aerogen
import AerodynamicsAirfoilNACA4Series as aerofoil4
import AerodynamicsUtil as aeroutil

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

epsilon = 1e-5

class Panel:
    def __init__(self):
        # Points should be ordered as:
        #
        # y  2----3
        # ^  |    |
        # |  |    |
        # |  1----4
        # |
        # +-----> x
        self.corners = np.zeros([4, 3])
        self.collacation = np.zeros(3)
        self.normal = np.zeros(3)
        self.xaxis = np.zeros(3)
        self.yaxis = np.zeros(3)
        self.indices = np.zeros(3, dtype=np.int32)
        self.influenceCompound = None
        self.influenceWing = None
        self.downwashWing = None
        self.influenceWake = None
        self.downwashWake = None
        self.downwash = 0
        self.velocity = np.zeros(3)
        self.strength = None
        self.wakeIndex = -1

    def setGeometry(self, corners, normal, shift=True):
        # copy corners and normal
        self.corners = corners
        self.normal = normal

        # construct x and y axis
        self.xaxis = self.corners[3] - self.corners[0]
        xaxisLength = np.linalg.norm(self.xaxis)

        if self.xaxis[1] != 0:
            raise ValueError("Expected y-component of constructed x-axis to be 0")

        self.yaxis = np.cross(self.normal, self.xaxis)
        self.xaxis /= xaxisLength
        self.yaxis /= np.linalg.norm(self.yaxis)

        # check if we should shift the panel along the x-axis (move the 1-2
        # vortex to the 0.25 panel position)
        if shift == True:
            for i in range(0, 4):
                self.corners[i] += self.xaxis * xaxisLength * 0.25

        # Calculate collacation point
        self.collacation = (self.corners[0] + self.corners[1] +
                            self.corners[2] + self.corners[3]) / 4.0

    def determineInfluence(self, wingPanels, wakePanels):
        self.influenceTotal = np.zeros(len(wingPanels))
        self.influenceWing = np.zeros(len(wingPanels))
        self.influenceWake = np.zeros(len(wakePanels))
        self.downwashWing = np.zeros(len(wingPanels))
        self.downwashWake = np.zeros(len(wakePanels))

        for i in range(0, len(wingPanels)):
            self.influenceTotal[i] = np.dot(wingPanels[i].determineVortexPanelInfluence(self.collacation), self.normal)
            self.influenceWing[i] = self.influenceTotal[i]
            self.downwashWing[i] = np.dot(wingPanels[i].determineVortexPanelDownwash(self.collacation), self.normal)

            if wingPanels[i].wakeIndex != -1:
                self.influenceTotal[i] += np.dot(wakePanels[wingPanels[i].wakeIndex].determineVortexPanelInfluence(self.collacation), self.normal)

        for i in range(0, len(wakePanels)):
            self.influenceWake[i] = np.dot(wakePanels[i].determineVortexPanelInfluence(self.collacation), self.normal)
            self.downwashWake[i] = np.dot(wakePanels[i].determineVortexPanelDownwash(self.collacation), self.normal)

    def determineVortexPanelInfluence(self, target):
        return (
            self.determineVortexLineInfluence(self.corners[0], self.corners[1], target) +
            self.determineVortexLineInfluence(self.corners[1], self.corners[2], target) +
            self.determineVortexLineInfluence(self.corners[2], self.corners[3], target) +
            self.determineVortexLineInfluence(self.corners[3], self.corners[0], target)
        )

    def determineVortexPanelDownwash(self, target):
        return (
            self.determineVortexLineInfluence(self.corners[1], self.corners[2], target) +
            self.determineVortexLineInfluence(self.corners[3], self.corners[0], target)
        )
    def determineVortexLineInfluence(self, pointA, pointB, target):
        ra = target - pointA
        rb = target - pointB
        rl = pointB - pointA

        #crossAB = np.cross(ra, rb)
        crossAB = np.asarray([
            (target[1] - pointA[1]) * (target[2] - pointB[2]) - (target[2] - pointA[2]) * (target[1] - pointB[1]),
            -(target[0] - pointA[0]) * (target[2] - pointB[2]) - (target[2] - pointA[2]) * (target[0] - pointB[0]),
            (target[0] - pointA[0]) * (target[1] - pointB[1]) - (target[1] - pointA[1]) * (target[0] - pointB[0])
        ])
        normA = np.linalg.norm(ra)
        normB = np.linalg.norm(rb)
        normAB = np.linalg.norm(crossAB)

        if normA < epsilon or normB < epsilon or normAB < epsilon:
            # Singularity
            return np.asarray([0, 0, 0])

        return crossAB * (1 / (4 * np.pi * normAB)) * (
            np.dot(rl, ra) / normA - np.dot(rl, rb) / normB
        )

class Analysis:
    def __init__(self, wakeEnd):
        self.wakeEnd = wakeEnd
        self.wings = []
        self.wingPanels = []
        self.wakePanels = []
        self.numSpan = []
        self.numChord = []

    def addWing(self, wing):
        wingPanels, wakePanels, numSpan, numChord = self.__generatePanels__(wing, len(self.wings))
        self.wings.append(wing)
        self.wingPanels.extend(wingPanels)
        self.wakePanels.extend(wakePanels)
        self.numSpan.append(numSpan)
        self.numChord.append(numChord)

    def analyze(self, Vinf):
        # Solve main matrix equation
        solutionMatrix = np.zeros([len(self.wingPanels), len(self.wingPanels)])
        downwashMatrix = np.zeros(solutionMatrix.shape)
        solutionVector = np.zeros(len(self.wingPanels))

        for i in range(0, len(self.wingPanels)):
            self.wingPanels[i].determineInfluence(self.wingPanels, self.wakePanels)

            for j in range(0, len(self.wingPanels)):
                solutionMatrix[i, j] = self.wingPanels[i].influenceWing[j]
                downwashMatrix[i, j] = self.wingPanels[i].downwashWing[j]

            solutionVector[i] = -np.dot(Vinf, self.wingPanels[i].normal)

        strengths = np.linalg.solve(solutionMatrix, solutionVector)
        downwash = np.matmul(downwashMatrix, strengths)

        for i in range(0, len(self.wingPanels)):
            self.wingPanels[i].strength = strengths[i]
            self.wingPanels[i].downwash = downwash[i]
            
            if self.wingPanels[i].wakeIndex != -1:
                self.wakePanels[self.wingPanels[i].wakeIndex].strength = strengths[i]
                
        # Calculate deltaLift
        self.deltaLift = np.zeros([self.numChord[0], self.numSpan[0]])
        for i in range(0, self.numSpan[0]):
            for j in range(0, self.numChord[0]):
                deltaY = 0
                curIndex = j * self.numSpan[0] + i
                if i == self.numSpan[0] - 1:
                    deltaY = self.wingPanels[curIndex].collacation[1] - \
                        self.wingPanels[curIndex - 1].collacation[1]
                else:
                    deltaY = self.wingPanels[curIndex + 1].collacation[1] - \
                        self.wingPanels[curIndex].collacation[1]

                if j == 0:
                    self.deltaLift[j, i] = self.wingPanels[curIndex].strength * deltaY
                else:
                    self.deltaLift[j, i] = (self.wingPanels[curIndex].strength -
                        self.wingPanels[(j - 1) * self.numSpan[0] + i].strength) * deltaY

    def __generatePanels__(self, wing, wingIndex):
        numPanelsSpan = wing.GetNumCamberPanelsY()
        numPanelsChord = wing.GetNumCamberPanelsX()

        wingPanels = [None] * numPanelsSpan * numPanelsChord
        wakePanels = [None] * numPanelsSpan

        # Handle all panels except for the last spanwise row
        for iSpan in range(0, numPanelsSpan):
            for iChord in range(0, numPanelsChord - 1):
                # Create panel
                newPanel = Panel()

                corners = np.zeros([4, 3])
                corners[0] = wing.GetCamberPanelBR(iChord, iSpan)
                corners[1] = wing.GetCamberPanelBL(iChord, iSpan)
                corners[2] = wing.GetCamberPanelTL(iChord, iSpan)
                corners[3] = wing.GetCamberPanelTR(iChord, iSpan)

                newPanel.setGeometry(corners, wing.GetCamberPanelNormal(iChord, iSpan))

                # Store new panel
                newPanel.indices = [wingIndex, iSpan, iChord]
                wingPanels[iChord * numPanelsSpan + iSpan] = newPanel

        # Process the last spanwise row and construct the wake at the same time
        for iSpan in range(0, numPanelsSpan):
            # Create wing panel
            newWingPanel = Panel()

            corners = np.zeros([4, 3])
            corners[0] = wing.GetCamberPanelBR(numPanelsChord - 1, iSpan)
            corners[1] = wing.GetCamberPanelBL(numPanelsChord - 1, iSpan)
            corners[2] = wing.GetCamberPanelTL(numPanelsChord - 1, iSpan)
            corners[3] = wing.GetCamberPanelTR(numPanelsChord - 1, iSpan)

            newWingPanel.setGeometry(corners, wing.GetCamberPanelNormal(numPanelsChord - 1, iSpan))

            # Store last wing panel
            newWingPanel.indices = [wingIndex, iSpan, numPanelsChord - 1]

            # Create wake panel by extending the last panel's trailing edge
            # to the specified limit
            newWakePanel = Panel()

            corners = np.zeros([4, 3])
            corners[0] = newWingPanel.corners[3]
            corners[1] = newWingPanel.corners[2]
            corners[2] = corners[0]
            corners[3] = corners[0]
            corners[2, 0] = self.wakeEnd
            corners[3, 0] = self.wakeEnd

            # Construct normal vector (as the x-direction changed)
            cross = np.cross(corners[2] - corners[0], corners[1] - corners[3])
            cross /= np.linalg.norm(cross)

            # Store new wake and wing panel
            newWakePanel.setGeometry(corners, cross, False)
            wakePanels[iSpan] = newWakePanel
            newWingPanel.wakeIndex = iSpan
            wingPanels[(numPanelsChord - 1) * numPanelsSpan + iSpan] = newWingPanel

        return wingPanels, wakePanels, numPanelsSpan, numPanelsChord

def __testAnalysisGeometry__():
    foil = aerofoil4.AirfoilNACA4Series('3412', aerogen.GeneratorDoubleChebyshev(0, 1, 20))
    wing = aerowing.Wing(foil, aerogen.GeneratorPiecewiseLinear([0.6, 2.0, 0.6], [8, 8]),
        aerogen.GeneratorLinear(0, 10, 16), 0.0, 0, [0, -5, 0])

    a = Analysis(20.0)
    a.addWing(wing)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Retrieve wing data
    wingPoints = np.zeros([len(a.wingPanels), 3])
    wingAxes = np.zeros([3, len(a.wingPanels), 3])

    for i in range(0, wingPoints.shape[0]):
        wingPoints[i] = a.wingPanels[i].collacation
        wingAxes[0, i] = a.wingPanels[i].xaxis
        wingAxes[1, i] = a.wingPanels[i].yaxis
        wingAxes[2, i] = a.wingPanels[i].normal

    ax.scatter(wingPoints[:,0], wingPoints[:,1], wingPoints[:,2], 'k')

    colors = ['r', 'g', 'b']
    for i in range(0, 3):
        ax.quiver(wingPoints[:,0], wingPoints[:,1], wingPoints[:,2],
                  wingAxes[i,:,0], wingAxes[i,:,1], wingAxes[i,:,2], length=0.1, color=colors[i])

    # Retrieve wake data
    wakePoints = np.zeros([len(a.wakePanels), 3])
    wakeAxes = np.zeros([3, len(a.wakePanels), 3])

    for i in range(0, wakePoints.shape[0]):
        wakePoints[i] = a.wakePanels[i].collacation
        wakeAxes[0, i] = a.wakePanels[i].xaxis
        wakeAxes[1, i] = a.wakePanels[i].yaxis
        wakeAxes[2, i] = a.wakePanels[i].normal

    ax.scatter(wakePoints[:, 0], wakePoints[:, 1], wakePoints[:, 2], 'k')
    for i in range(0, 3):
        ax.quiver(wakePoints[:,0], wakePoints[:,1], wakePoints[:,2],
                  wakeAxes[i,:,0], wakeAxes[i,:,1], wakeAxes[i,:,2], length=0.1, color=colors[i])

    aeroutil.set3DAxesEqual(ax, wingPoints[:,0], wingPoints[:,1], wingPoints[:,2])

def __testAnalysisParameters__():
    foil = aerofoil4.AirfoilNACA4Series('9436', aerogen.GeneratorDoubleChebyshev(0, 1, 25))
    wing = aerowing.Wing(foil, aerogen.GeneratorPiecewiseLinear([1.5, 3.5, 1.5], [15, 16]),
        aerogen.GeneratorLinear(0, 15, 31), 0.0, 0, [0, -7.5, 0])

    a = Analysis(20.0)
    a.addWing(wing)
    a.analyze([100, 0, 30])

    # Retrieve the triangulation
    coords = wing.camberCoordinates
    triang = wing.camberTriangulation
    values = np.zeros([len(triang)])
    colors = np.zeros([len(triang), 4])

    cmap = plt.get_cmap("jet")

    totalLift = np.zeros(a.numSpan[0])
    for i in range(0, a.numChord[0]):
        for j in range(0, a.numSpan[0]):
            indexWing = j * a.numChord[0] + i
            values[2 * indexWing] = a.deltaLift[i, j]
            values[2 * indexWing + 1] = a.deltaLift[i, j]
            #values[2 * indexWing] = a.wingPanels[i * a.numSpan[0] + j].downwash
            #values[2 * indexWing + 1] = values[2 * indexWing]
            if i == 1:
                values[2 * (indexWing - 1)] = values[2 * indexWing]
                values[2 * (indexWing - 1) + 1] = values[2 * indexWing]
                totalLift[j] = a.deltaLift[i, j]
            #values[2 * indexWing] = a.wingPanels[i * a.numSpan[0] + j].strength
            #values[2 * indexWing + 1] = values[2 * indexWing]
            totalLift[j] += a.deltaLift[i, j]
            
    print(totalLift)

    minValue = min(values)
    maxValue = max(values)
    
    print('min =', minValue)
    print('max =', maxValue)

    for i in range(0, len(values)):
        colors[i] = cmap((values[i] - minValue) / (maxValue - minValue))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    aeroutil.plot3DColors(ax, coords[:, 0], coords[:, 1], coords[:, 2], triang, colors)
    aeroutil.set3DAxesEqual(ax, coords[:, 0], coords[:, 1], coords[:, 2])

def __testAnalysis__():
    foil = aerofoil4.AirfoilNACA4Series('8436', aerogen.GeneratorDoubleChebyshev(0, 1, 20))
    wing = aerowing.Wing(foil, aerogen.GeneratorPiecewiseLinear([0.6, 2.0, 0.6], [8, 8]),
        aerogen.GeneratorLinear(0, 10, 16), 0.0, 0, [0, -5, 0])

    a = Analysis(20.0)
    a.addWing(wing)
    a.analyze([100, 0, 20])

    # Retrieve wing data
    wingPoints = np.zeros([len(a.wingPanels), 3])
    wingVelocity = np.zeros(wingPoints.shape)

    for i in range(0, len(a.wingPanels)):
        wingPoints[i] = a.wingPanels[i].collacation
        wingVelocity[i] = a.wingPanels[i].velocity

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(wingPoints[:, 0], wingPoints[:, 1], wingPoints[:, 2], 'k')
    ax.quiver(wingPoints[:, 0], wingPoints[:, 1], wingPoints[:, 2],
              wingVelocity[:, 0], wingVelocity[:, 1], wingVelocity[:, 2], length=0.25)

    print(wingVelocity)

    aeroutil.set3DAxesEqual(ax, wingPoints[:, 0], wingPoints[:, 1], wingPoints[:, 2])

__testAnalysisParameters__()
