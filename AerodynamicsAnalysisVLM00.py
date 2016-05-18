# -*- coding: utf-8 -*-
"""
Created on Tue May 17 17:02:34 2016

@author: MaxHenger
"""

import AerodynamicsWing as aerowing
import AerodynamicsAirfoilNACA4Series as aerofoil4
import AerodynamicsGenerator as aerogen
import AerodynamicsUtil as aeroutil

import numpy as np
from enum import Enum
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class AnalysisIndex(Enum):
    CollacationPoint = 0
    Corners = 1
    Area = 2
    NormalVector = 3
    XVector = 4
    YVector = 5
    SourceStrength = 6
    DoubletStrength = 7
    InfluenceWingDoublet = 8
    InfluenceWingSource = 9
    InfluenceWakeDoublet = 10
    LocalVelocityWingDoublet = 11
    LocalVelocityWingSource = 12
    LocalVelocityWakeDoublet = 13
    Velocity = 14
    IsUpper = 15
    WakePanel = 16
    Indices = 17
    Total = 18

class AnalysisVLM00:
    def __init__(self, numWake, wakeEnd):
        self.wings = []
        self.wingsChanged = True
        self.numWake = numWake
        self.wakeEnd = wakeEnd
        self.wingPanels = []
        self.wakePanels = []

    def AddWing(self, wing):
        # Check if the wing is the expected AerodynamicsWing.Wing object
        if not isinstance(wing, aerowing.Wing):
            raise ValueError("Expected variable 'wing' to be an instance of the " +
                "AerodynamicsWing.Wing class")

        if wing.GetNumUpperPanelsY() != wing.GetNumLowerPanelsY():
            raise ValueError("Expected wing to contain the same number of panels " +
                "along the y-axis on the upper side as on the lower side")

        self.wings.append(wing)

    def Analyze(self, Vinf):
        # Update influence parameters if necessary
        self.Vinf = Vinf

        if self.wingsChanged == True:
            self.__updateWingGeometry__()
            self.__updateInfluenceParameters__()

        # Assume the user knows what he/she is doing: calculate singularity
        # strengths
        self.__updateSingularityStrengths__(Vinf)

        # Calculate the velocity offset at each panel
        for iWingConsidered in range(0, len(self.wingPanels)):
            curPanel = self.wingPanels[iWingConsidered]
            curPanel[AnalysisIndex.Velocity.value] = Vinf

            for iWing in range(0, len(self.wingPanels)):
                curPanel[AnalysisIndex.Velocity.value] += \
                    self.wingPanels[iWing][AnalysisIndex.SourceStrength.value] * \
                    curPanel[AnalysisIndex.LocalVelocityWingSource.value][iWing]

                curPanel[AnalysisIndex.Velocity.value] += \
                    self.wingPanels[iWing][AnalysisIndex.DoubletStrength.value] * \
                    curPanel[AnalysisIndex.LocalVelocityWingDoublet.value][iWing]

            for iWake in range(0, len(self.wakePanels)):
                curPanel[AnalysisIndex.Velocity.value] += \
                    self.wakePanels[iWake][AnalysisIndex.DoubletStrength.value] * \
                    curPanel[AnalysisIndex.LocalVelocityWakeDoublet.value][iWake]

    def GetParameters(self, iWing):
        curWing = self.wings[iWing]
        numChordLower = curWing.GetNumLowerPanelsX()
        numChordUpper = curWing.GetNumUpperPanelsX()

        paramsLower = [None] * numChordLower * curWing.GetNumLowerPanelsY()
        paramsUpper = [None] * numChordUpper * curWing.GetNumUpperPanelsY()

        for i in range(0, len(self.wingPanels)):
            curPanel = self.wingPanels[i]
            curIndices = curPanel[AnalysisIndex.Indices.value]

            if curIndices[0] == iWing:
                if curIndices[3] == True:
                    paramsUpper[curIndices[2] * numChordUpper + curIndices[1]] = curPanel
                else:
                    paramsLower[curIndices[2] * numChordLower + curIndices[1]] = curPanel

        return paramsUpper, paramsLower

    def __updateWingGeometry__(self):
        self.wingPanels = []
        self.wakePanels = []

        for iWing in range(0, len(self.wings)):
            localWingPanels, localWakePanels = self.__generateWingGeometry__(iWing, self.wings[iWing], len(self.wakePanels))
            self.wingPanels.extend(localWingPanels)
            self.wakePanels.extend(localWakePanels)

    def __constructXYAxes__(self, center, utilPoint1, utilPoint2, normal):
        localX = 0.5 * (utilPoint1 + utilPoint2) - center
        localY = np.cross(normal, localX)
        return localX / np.linalg.norm(localX), localY / np.linalg.norm(localY)

    def __generateWingGeometry__(self, wingIndex, wing, wakeOffset):
        # Convert wing data into a set of upper surface-, lower surface- and
        # wake panels
        wakePanels = []
        wingPanels = []

        numX = wing.GetNumLowerPanelsX()
        numY = wing.GetNumLowerPanelsY()

        # Move along span, do each row of upper and lower panels while
        # generating the wake. It is assumed that the airfoil closes neatly
        # such that a wake only has to be attached to one wake panel.
        # Keep track of the created wake panels such that the upper panels can
        # be later linked to them as well
        baseWake = []

        for iY in range(0, numY):
            for iX in range(0, numX - 1):
                # Generate new lower wing panel
                newPanel = [None] * AnalysisIndex.Total.value

                center = wing.GetLowerPanelCollacationPoint(iX, iY)
                points = wing.GetLowerPanelPoints(iX, iY)
                normal = wing.GetLowerPanelNormal(iX, iY)

                newPanel[AnalysisIndex.CollacationPoint.value] = center
                newPanel[AnalysisIndex.Corners.value] = points
                newPanel[AnalysisIndex.Area.value] = wing.GetLowerPanelArea(iX, iY)
                newPanel[AnalysisIndex.NormalVector.value] = normal
                newPanel[AnalysisIndex.IsUpper.value] = False
                newPanel[AnalysisIndex.WakePanel.value] = -1
                newPanel[AnalysisIndex.XVector.value], newPanel[AnalysisIndex.YVector.value] = \
                    self.__constructXYAxes__(center, points[0], points[1], normal)
                newPanel[AnalysisIndex.Indices.value] = [wingIndex, iX, iY, False]

                wingPanels.append(newPanel)

            # Add the last panel, which is connected to the wake
            newPanel = [None] * AnalysisIndex.Total.value

            center = wing.GetLowerPanelCollacationPoint(numX - 1, iY)
            points = wing.GetLowerPanelPoints(numX - 1, iY)
            normal = wing.GetLowerPanelNormal(numX - 1, iY)

            newPanel[AnalysisIndex.CollacationPoint.value] = center
            newPanel[AnalysisIndex.Corners.value] = points
            newPanel[AnalysisIndex.Area.value] = wing.GetLowerPanelArea(numX - 1, iY)
            newPanel[AnalysisIndex.NormalVector.value] = normal
            newPanel[AnalysisIndex.IsUpper.value] = False
            newPanel[AnalysisIndex.WakePanel.value] = wakeOffset + len(wakePanels)
            newPanel[AnalysisIndex.XVector.value], newPanel[AnalysisIndex.YVector.value] = \
                self.__constructXYAxes__(center, points[0], points[1], normal)
            newPanel[AnalysisIndex.Indices.value] = [wingIndex, numX - 1, iY, False]

            wingPanels.append(newPanel)
            baseWake.append(wakeOffset + len(wakePanels))

            # Generate variables used to generate wake panels
            tempCoords = wing.GetLowerPanelPoints(iX, iY)
            pointInner = tempCoords[0]
            pointOuter = tempCoords[1]
            interpolateInner = np.linspace(pointInner[0], self.wakeEnd, self.numWake + 1)
            interpolateOuter = np.linspace(pointOuter[0], self.wakeEnd, self.numWake + 1)

            # Create all wake panels
            for iW in range(0, self.numWake):
                newPanel = [None] * AnalysisIndex.Total.value

                # Construct corners and the collacation point
                newPanel[AnalysisIndex.Corners.value] = np.asarray([
                    [interpolateInner[iW], pointInner[1], pointInner[2]],
                    [interpolateInner[iW + 1], pointInner[1], pointInner[2]],
                    [interpolateOuter[iW + 1], pointOuter[1], pointOuter[2]],
                    [interpolateOuter[iW], pointOuter[1], pointOuter[2]]
                ])
                newPanel[AnalysisIndex.CollacationPoint.value] = (
                    newPanel[AnalysisIndex.Corners.value][0] +
                    newPanel[AnalysisIndex.Corners.value][1] +
                    newPanel[AnalysisIndex.Corners.value][2] +
                    newPanel[AnalysisIndex.Corners.value][3]
                ) / 4.0

                # Generate normal vector and area
                toTL = newPanel[AnalysisIndex.Corners.value][2] - newPanel[AnalysisIndex.Corners.value][3]
                toTR = newPanel[AnalysisIndex.Corners.value][1] - newPanel[AnalysisIndex.Corners.value][3]
                toBR = newPanel[AnalysisIndex.Corners.value][0] - newPanel[AnalysisIndex.Corners.value][3]

                crossTLAndTR = np.cross(toTR, toTL)
                crossTRAndBR = np.cross(toBR, toTR)

                normal = crossTLAndTR + crossTRAndBR
                normal /= np.linalg.norm(normal)

                newPanel[AnalysisIndex.NormalVector.value] = normal
                newPanel[AnalysisIndex.Area.value] = 0.5 * (np.linalg.norm(crossTLAndTR) + np.linalg.norm(crossTRAndBR))

                # Generate local X and Y vector
                newPanel[AnalysisIndex.XVector.value], newPanel[AnalysisIndex.YVector.value] = \
                    self.__constructXYAxes__(newPanel[AnalysisIndex.CollacationPoint.value],
                                             newPanel[AnalysisIndex.Corners.value][0],
                                             newPanel[AnalysisIndex.Corners.value][1],
                                             newPanel[AnalysisIndex.NormalVector.value])

                wakePanels.append(newPanel)

        numX = wing.GetNumUpperPanelsX()
        numY = wing.GetNumUpperPanelsY()

        for iY in range(0, numY):
            for iX in range(0, numX - 1):
                # Generate new upper wing panel
                newPanel = [None] * AnalysisIndex.Total.value

                center = wing.GetUpperPanelCollacationPoint(iX, iY)
                points = wing.GetUpperPanelPoints(iX, iY)
                normal = wing.GetUpperPanelNormal(iX, iY)

                newPanel[AnalysisIndex.CollacationPoint.value] = center
                newPanel[AnalysisIndex.Corners.value] = points
                newPanel[AnalysisIndex.Area.value] = wing.GetUpperPanelArea(iX, iY)
                newPanel[AnalysisIndex.NormalVector.value] = normal
                newPanel[AnalysisIndex.IsUpper.value] = True
                newPanel[AnalysisIndex.WakePanel.value] = -1
                newPanel[AnalysisIndex.XVector.value], newPanel[AnalysisIndex.YVector.value] = \
                    self.__constructXYAxes__(center, points[0], points[1], normal)
                newPanel[AnalysisIndex.Indices.value] = [wingIndex, iX, iY, True]

                wingPanels.append(newPanel)

            # Generate last upper panel connected to the wake
            newPanel = [None] * AnalysisIndex.Total.value

            center = wing.GetUpperPanelCollacationPoint(numX - 1, iY)
            points = wing.GetUpperPanelPoints(numX - 1, iY)
            normal = wing.GetUpperPanelNormal(numX - 1, iY)

            newPanel[AnalysisIndex.CollacationPoint.value] = center
            newPanel[AnalysisIndex.Corners.value] = points
            newPanel[AnalysisIndex.NormalVector.value] = normal
            newPanel[AnalysisIndex.IsUpper.value] = True
            newPanel[AnalysisIndex.WakePanel.value] = baseWake[iY]
            newPanel[AnalysisIndex.XVector.value], newPanel[AnalysisIndex.YVector.value] = \
                self.__constructXYAxes__(center, points[0], points[1], normal)
            newPanel[AnalysisIndex.Indices.value] = [wingIndex, numX - 1, iY, True]

            wingPanels.append(newPanel)

        return wingPanels, wakePanels

    def __updateInfluenceParameters__(self):
        # For each wing panel calculate the influences of all other panels
        for iWingConsidered in range(0, len(self.wingPanels)):
            # Determine influence of other wing panel source and doublet
            # singularity functions
            curPanel = self.wingPanels[iWingConsidered]
            curPanel[AnalysisIndex.InfluenceWingSource.value] = np.zeros(len(self.wingPanels))
            curPanel[AnalysisIndex.LocalVelocityWingSource.value] = np.zeros([len(self.wingPanels), 3])
            curPanel[AnalysisIndex.InfluenceWingDoublet.value] = np.zeros(len(self.wingPanels))
            curPanel[AnalysisIndex.LocalVelocityWingDoublet.value] = np.zeros([len(self.wingPanels), 3])

            for iWingAffecting in range(0, len(self.wingPanels)):
                curPanel[AnalysisIndex.InfluenceWingSource.value][iWingAffecting], \
                curPanel[AnalysisIndex.LocalVelocityWingSource.value][iWingAffecting] = \
                    self.__calculateSourcePanelInfluence__(
                        curPanel[AnalysisIndex.CollacationPoint.value],
                        curPanel[AnalysisIndex.Corners.value],
                        curPanel[AnalysisIndex.NormalVector.value],
                        curPanel[AnalysisIndex.XVector.value],
                        curPanel[AnalysisIndex.YVector.value],
                        self.wingPanels[iWingAffecting][AnalysisIndex.CollacationPoint.value]
                    )

                curPanel[AnalysisIndex.InfluenceWingDoublet.value][iWingAffecting], \
                curPanel[AnalysisIndex.LocalVelocityWingDoublet.value][iWingAffecting]= \
                    self.__calculateDoubletPanelInfluence__(
                        curPanel[AnalysisIndex.CollacationPoint.value],
                        curPanel[AnalysisIndex.Corners.value],
                        curPanel[AnalysisIndex.NormalVector.value],
                        curPanel[AnalysisIndex.XVector.value],
                        curPanel[AnalysisIndex.YVector.value],
                        self.wingPanels[iWingAffecting][AnalysisIndex.CollacationPoint.value]
                    )

            # Determine influence of wake panels
            curPanel[AnalysisIndex.InfluenceWakeDoublet.value] = np.zeros(len(self.wakePanels))
            curPanel[AnalysisIndex.LocalVelocityWakeDoublet.value] = np.zeros([len(self.wakePanels), 3])

            for iWakeAffecting in range(0, len(self.wakePanels)):
                curPanel[AnalysisIndex.InfluenceWakeDoublet.value][iWakeAffecting], \
                curPanel[AnalysisIndex.LocalVelocityWakeDoublet.value][iWakeAffecting] = \
                    self.__calculateDoubletPanelInfluence__(
                        curPanel[AnalysisIndex.CollacationPoint.value],
                        curPanel[AnalysisIndex.Corners.value],
                        curPanel[AnalysisIndex.NormalVector.value],
                        curPanel[AnalysisIndex.XVector.value],
                        curPanel[AnalysisIndex.YVector.value],
                        self.wakePanels[iWakeAffecting][AnalysisIndex.CollacationPoint.value]
                    )

    def __updateSingularityStrengths__(self, Vinf):
        # Update all wing panel source strengths
        for iWing in range(0, len(self.wingPanels)):
            self.wingPanels[iWing][AnalysisIndex.SourceStrength.value] = \
                - np.dot(self.wingPanels[iWing][AnalysisIndex.NormalVector.value], Vinf)

        # Generate the matrix equation that will be solved
        matrix = np.zeros([len(self.wingPanels), len(self.wingPanels)])
        vector = np.zeros(len(self.wingPanels))
        wakeLookup = {}

        for iWingConsidered in range(0, len(self.wingPanels)):
            curWing = self.wingPanels[iWingConsidered]

            for iWingAffecting in range(0, len(self.wingPanels)):
                curAffecting = self.wingPanels[iWingAffecting]

                vector[iWingConsidered] -= \
                    curWing[AnalysisIndex.InfluenceWingSource.value][iWingAffecting] * \
                    curAffecting[AnalysisIndex.SourceStrength.value]

                if curAffecting[AnalysisIndex.WakePanel.value] != -1:
                    # A wake is connected to this panel
                    if self.wingPanels[iWing][AnalysisIndex.IsUpper.value] == True:
                        matrix[iWingConsidered, iWingAffecting] += curWing[AnalysisIndex.InfluenceWingDoublet.value][iWingAffecting] + \
                            curWing[AnalysisIndex.InfluenceWakeDoublet.value][curAffecting[AnalysisIndex.WakePanel.value]]
                    else:
                        matrix[iWingConsidered, iWingAffecting] += curWing[AnalysisIndex.InfluenceWingDoublet.value][iWingAffecting] - \
                            curWing[AnalysisIndex.InfluenceWakeDoublet.value][curAffecting[AnalysisIndex.WakePanel.value]]
                else:
                    matrix[iWingConsidered, iWingAffecting] += curWing[AnalysisIndex.InfluenceWingDoublet.value][iWingAffecting]

            curWakeIndex = curWing[AnalysisIndex.WakePanel.value]
            if curWakeIndex != -1:
                # Update the wake lookup
                if curWakeIndex in wakeLookup.keys():
                    if curWing[AnalysisIndex.IsUpper.value] == True:
                        wakeLookup[curWakeIndex][0] = iWingConsidered
                    else:
                        wakeLookup[curWakeIndex][1] = iWingConsidered
                else:
                    if curWing[AnalysisIndex.IsUpper.value] == True:
                        wakeLookup[curWakeIndex] = [iWingConsidered, -1]
                    else:
                        wakeLookup[curWakeIndex] = [-1, iWingConsidered]

        # Solve the matrix equation to obtain the wing doublet strengths
        solution = np.linalg.solve(matrix, vector)

        for iWingConsidered in range(0, len(self.wingPanels)):
            curWing = self.wingPanels[iWingConsidered]
            curWing[AnalysisIndex.DoubletStrength.value] = solution[iWingConsidered]

        for wakeBaseIndex, associatedPanels in wakeLookup.items():
            doubletStrength = self.wingPanels[associatedPanels[0]][AnalysisIndex.DoubletStrength.value] - \
                self.wingPanels[associatedPanels[1]][AnalysisIndex.DoubletStrength.value]

            for iWake in range(wakeBaseIndex, wakeBaseIndex + self.numWake):
                self.wakePanels[iWake][AnalysisIndex.DoubletStrength.value] = doubletStrength

    def __calculateAuxilliary__(self, panelCenter, panelCoords, panelNormal, panelX, panelY, targetPoint):
        # Note that this algorithm is taken from "Low-Speed Aerodynamics: From
        # Wing Theory to Panel Methods" by Joseph Katz and Allen Plotkin

        # Move all coordinates to the local reference frame
        panelCoordinates = np.asarray([
            panelCoords[0] - panelCenter,
            panelCoords[1] - panelCenter,
            panelCoords[2] - panelCenter,
            panelCoords[3] - panelCenter
        ])
        point = targetPoint - panelCenter

        # Project the point onto the axes to get a relative vector with respect
        # to the panel's reference frame
        localPoint = [
            np.dot(point, panelX),
            np.dot(point, panelY),
            np.dot(point, panelNormal)
        ]

        # Calculate auxilliary variables
        r = np.zeros(4)
        e = np.zeros(4)
        h = np.zeros(4)

        for i in range(0, 4):
            r[i] = np.sqrt(
                np.power(localPoint[0] - panelCoordinates[i, 0], 2.0) +
                np.power(localPoint[1] - panelCoordinates[i, 1], 2.0) +
                np.power(localPoint[2], 2.0)
            )
            e[i] = np.power(localPoint[0] - panelCoordinates[i, 0], 2.0) + \
                   np.power(localPoint[2], 2.0)
            h[i] = (localPoint[0] - panelCoordinates[i, 0]) * \
                   (localPoint[1] - panelCoordinates[i, 1])

        m12 = (panelCoordinates[1, 1] - panelCoordinates[0, 1]) / \
            (panelCoordinates[1, 0] - panelCoordinates[0, 0])
        m23 = (panelCoordinates[2, 1] - panelCoordinates[1, 1]) / \
            (panelCoordinates[2, 0] - panelCoordinates[1, 0])
        m34 = (panelCoordinates[3, 1] - panelCoordinates[2, 1]) / \
            (panelCoordinates[3, 0] - panelCoordinates[2, 0])
        m41 = (panelCoordinates[0, 1] - panelCoordinates[3, 1]) / \
            (panelCoordinates[0, 0] - panelCoordinates[3, 0])

        d12 = np.sqrt(
            np.power(panelCoordinates[1, 0] - panelCoordinates[0, 0], 2.0) +
            np.power(panelCoordinates[1, 1] - panelCoordinates[0, 1], 2.0)
        )
        d23 = np.sqrt(
            np.power(panelCoordinates[2, 0] - panelCoordinates[1, 0], 2.0) +
            np.power(panelCoordinates[2, 1] - panelCoordinates[1, 1], 2.0)
        )
        d34 = np.sqrt(
            np.power(panelCoordinates[3, 0] - panelCoordinates[2, 0], 2.0) +
            np.power(panelCoordinates[3, 1] - panelCoordinates[2, 1], 2.0)
        )
        d41 = np.sqrt(
            np.power(panelCoordinates[0, 0] - panelCoordinates[3, 0], 2.0) +
            np.power(panelCoordinates[0, 1] - panelCoordinates[3, 1], 2.0)
        )

        return (panelCoordinates, localPoint,
                r, e, h,
                m12, m23, m34, m41,
                d12, d23, d34, d41)

    def __calculateSourcePanelInfluence__(self, panelCenter, panelCoordinates,
                                          panelNormal, panelX, panelY, point):
        # Note that this algorithm is taken from "Low-Speed Aerodynamics: From
        # Wing Theory to Panel Methods" by Joseph Katz and Allen Plotkin

        # Retrieve auxilliary variables
        (panelCoordinates, localPoint, r, e, h,
         m12, m23, m34, m41, d12, d23, d34, d41) = \
            self.__calculateAuxilliary__(panelCenter, panelCoordinates, panelNormal,
                                         panelX, panelY, point)

        if r[0] < 1e-5 or r[1] < 1e-5:
            return (0, [0, 0, 0])

        # Calculate velocity components for a later transformation
        u = ( 1.0 / (4.0 * np.pi) * (
            (panelCoordinates[1, 1] - panelCoordinates[0, 1]) / d12 *
            np.log( (r[0] + r[1] - d12) / (r[0] + r[1] + d12) ) +
            (panelCoordinates[2, 1] - panelCoordinates[1, 1]) / d23 *
            np.log( (r[1] + r[2] - d23) / (r[1] + r[2] + d23) ) +
            (panelCoordinates[3, 1] - panelCoordinates[2, 1]) / d34 *
            np.log( (r[2] + r[3] - d34) / (r[2] + r[3] + d34) ) +
            (panelCoordinates[0, 1] - panelCoordinates[3, 1]) / d41 *
            np.log( (r[3] + r[0] - d41) / (r[3] + r[0] + d41) )
        ))
        v = ( 1.0 / (4.0 * np.pi) * (
            (panelCoordinates[0, 0] - panelCoordinates[1, 0]) / d12 *
            np.log( (r[0] + r[1] - d12) / (r[0] + r[1] + d12) ) +
            (panelCoordinates[1, 0] - panelCoordinates[2, 0]) / d23 *
            np.log( (r[1] + r[2] - d23) / (r[1] + r[2] + d23) ) +
            (panelCoordinates[2, 0] - panelCoordinates[3, 0]) / d34 *
            np.log( (r[2] + r[3] - d34) / (r[2] + r[3] + d34) ) +
            (panelCoordinates[3, 0] - panelCoordinates[0, 0]) / d41 *
            np.log( (r[3] + r[0] - d41) / (r[3] + r[0] + d41) )
        ))
        w = ( 1.0 / (4.0 * np.pi) * (
            np.arctan2(m12 * e[0] - h[0], localPoint[2] * r[0]) -
            np.arctan2(m12 * e[1] - h[1], localPoint[2] * r[1]) +
            np.arctan2(m23 * e[1] - h[1], localPoint[2] * r[1]) -
            np.arctan2(m23 * e[2] - h[2], localPoint[2] * r[2]) +
            np.arctan2(m34 * e[2] - h[2], localPoint[2] * r[2]) -
            np.arctan2(m34 * e[3] - h[3], localPoint[2] * r[3]) +
            np.arctan2(m41 * e[3] - h[3], localPoint[2] * r[3]) -
            np.arctan2(m41 * e[0] - h[0], localPoint[2] * r[0])
        ))

        # Calculate total influence
        return \
        (   # Influence parameters
            -1.0 / (4.0 * np.pi) * (
                (
                    (
                        (localPoint[0] - panelCoordinates[0, 0]) *
                        (panelCoordinates[1, 1] - panelCoordinates[0, 1]) -
                        (localPoint[1] - panelCoordinates[0, 1]) *
                        (panelCoordinates[1, 0] - panelCoordinates[0, 0])
                    ) / d12 * np.log( (r[0] + r[1] + d12) / (r[0] + r[1] - d12) ) +
                    (
                        (localPoint[0] - panelCoordinates[1, 0]) *
                        (panelCoordinates[2, 1] - panelCoordinates[1, 1]) -
                        (localPoint[1] - panelCoordinates[1, 1]) *
                        (panelCoordinates[2, 0] - panelCoordinates[1, 0])
                    ) / d23 * np.log( (r[1] + r[2] + d23) / (r[1] + r[2] - d23) ) +
                    (
                        (localPoint[0] - panelCoordinates[2, 0]) *
                        (panelCoordinates[3, 1] - panelCoordinates[2, 1]) -
                        (localPoint[1] - panelCoordinates[2, 1]) *
                        (panelCoordinates[3, 0] - panelCoordinates[2, 0])
                    ) / d34 * np.log( (r[2] + r[3] + d34) / (r[2] + r[3] - d34) ) +
                    (
                        (localPoint[0] - panelCoordinates[3, 0]) *
                        (panelCoordinates[0, 1] - panelCoordinates[3, 1]) -
                        (localPoint[1] - panelCoordinates[3, 1]) *
                        (panelCoordinates[0, 0] - panelCoordinates[3, 0])
                    ) / d41 * np.log( (r[3] + r[0] + d41) / (r[3] + r[0] - d41) )
                ) +
                abs(localPoint[2]) * (
                    np.arctan2(m12 * e[0] - h[0], localPoint[2] * r[0]) -
                    np.arctan2(m12 * e[1] - h[1], localPoint[2] * r[1]) +
                    np.arctan2(m23 * e[1] - h[1], localPoint[2] * r[1]) -
                    np.arctan2(m23 * e[2] - h[2], localPoint[2] * r[2]) +
                    np.arctan2(m34 * e[2] - h[2], localPoint[2] * r[2]) -
                    np.arctan2(m34 * e[3] - h[3], localPoint[2] * r[3]) +
                    np.arctan2(m41 * e[3] - h[3], localPoint[2] * r[3]) -
                    np.arctan2(m41 * e[0] - h[0], localPoint[2] * r[0])
                )
            ),
            u * panelX + v * panelY + w * panelNormal
        )

    def __calculateDoubletPanelInfluence__(self, panelCenter, panelCoordinates,
                                           panelNormal, panelX, panelY, point):
        # Note that this algorithm is taken from "Low-Speed Aerodynamics: From
        # Wing Theory to Panel Methods" by Joseph Katz and Allen Plotkin

        # Retrieve auxilliary variables
        (_, localPoint, r, e, h, m12, m23, m34, m41, d12, d23, d34, d41) = \
            self.__calculateAuxilliary__(panelCenter, panelCoordinates, panelNormal,
                                         panelX, panelY, point)

        if r[0] < 1e-5 or r[1] < 1e-5:
            return (0, [0, 0, 0])

        # Calculate often-used terms
        term1 = (r[0] * r[1]) / (r[0] * r[1] * (r[0] * r[1] - (
            (localPoint[0] - panelCoordinates[0, 0]) *
            (localPoint[0] - panelCoordinates[1, 0]) +
            (localPoint[1] - panelCoordinates[0, 1]) *
            (localPoint[1] - panelCoordinates[1, 1]) +
            np.power(localPoint[2], 2.0)
        )))
        term2 = (r[1] + r[2]) / (r[1] * r[2] * (r[1] * r[2] - (
            (localPoint[0] - panelCoordinates[1, 0]) *
            (localPoint[0] - panelCoordinates[2, 0]) +
            (localPoint[1] - panelCoordinates[1, 1]) *
            (localPoint[1] - panelCoordinates[2, 1]) +
            np.power(localPoint[2], 2.0)
        )))
        term3 = (r[2] + r[3]) / (r[2] * r[3] * (r[2] * r[3] - (
            (localPoint[0] - panelCoordinates[2, 0]) *
            (localPoint[0] - panelCoordinates[3, 0]) +
            (localPoint[1] - panelCoordinates[2, 1]) *
            (localPoint[1] - panelCoordinates[3, 1]) +
            np.power(localPoint[2], 2.0)
        )))
        term4 = (r[3] + r[0]) / (r[3] * r[0] * (r[3] * r[0] - (
            (localPoint[0] - panelCoordinates[3, 0]) *
            (localPoint[0] - panelCoordinates[0, 0]) +
            (localPoint[1] - panelCoordinates[3, 1]) *
            (localPoint[1] - panelCoordinates[0, 1]) +
            np.power(localPoint[2], 2.0)
        )))

        # Precalculate velocity components for a transformation
        u = 1.0 / (4.0 * np.pi) * (
            localPoint[2] * (panelCoordinates[0, 1] - panelCoordinates[1, 1]) * term1 +
            localPoint[2] * (panelCoordinates[1, 1] - panelCoordinates[2, 1]) * term2 +
            localPoint[2] * (panelCoordinates[2, 1] - panelCoordinates[3, 1]) * term3 +
            localPoint[2] * (panelCoordinates[3, 1] - panelCoordinates[0, 1]) * term4
        )
        v = 1.0 / (4.0 * np.pi) * (
            localPoint[2] * (panelCoordinates[1, 0] - panelCoordinates[0, 0]) * term1 +
            localPoint[2] * (panelCoordinates[2, 0] - panelCoordinates[1, 0]) * term2 +
            localPoint[2] * (panelCoordinates[3, 0] - panelCoordinates[2, 0]) * term3 +
            localPoint[2] * (panelCoordinates[0, 0] - panelCoordinates[3, 0]) * term4
        ),
        w = 1.0 / (4.0 * np.pi) * (
            (
                (localPoint[0] - panelCoordinates[1, 0]) *
                (localPoint[1] - panelCoordinates[0, 1]) -
                (localPoint[0] - panelCoordinates[0, 0]) *
                (localPoint[1] - panelCoordinates[1, 1])
            ) * term1 +
            (
                (localPoint[0] - panelCoordinates[2, 0]) *
                (localPoint[1] - panelCoordinates[1, 1]) -
                (localPoint[0] - panelCoordinates[1, 0]) *
                (localPoint[1] - panelCoordinates[2, 1])
            ) * term2 +
            (
                (localPoint[0] - panelCoordinates[3, 0]) *
                (localPoint[1] - panelCoordinates[2, 1]) -
                (localPoint[0] - panelCoordinates[2, 0]) *
                (localPoint[1] - panelCoordinates[3, 1])
            ) * term3 +
            (
                (localPoint[0] - panelCoordinates[0, 0]) *
                (localPoint[1] - panelCoordinates[3, 1]) -
                (localPoint[0] - panelCoordinates[3, 0]) *
                (localPoint[1] - panelCoordinates[0, 1])
            ) * term4
        )

        # Calculate total influence parameters
        return \
        (   # Influence parameters
            1.0 / (4.0 * np.pi) * (
                np.arctan2(m12 * e[0] * h[0], localPoint[2] * r[0]) -
                np.arctan2(m12 * e[1] - h[1], localPoint[2] * r[1]) +
                np.arctan2(m23 * e[1] - h[1], localPoint[2] * r[1]) -
                np.arctan2(m23 * e[2] - h[2], localPoint[2] * r[2]) +
                np.arctan2(m34 * e[2] - h[2], localPoint[2] * r[2]) -
                np.arctan2(m34 * e[3] - h[3], localPoint[2] * r[3]) +
                np.arctan2(m41 * e[3] - h[3], localPoint[2] * r[3]) -
                np.arctan2(m41 * e[0] - h[0], localPoint[2] * r[0])
            ),
            u * panelX + v * panelY + w * panelNormal
        )

def __testVLM00Geometry__():
    foil = aerofoil4.AirfoilNACA4Series('2436',
        aerogen.GeneratorSingleChebyshev(0.0, 1.0, 12))

    wing = aerowing.Wing(foil, aerogen.GeneratorLinear(2.0, 0.6, 10),
        aerogen.GeneratorLinear(0, 5, 10), np.pi / 12.0, 0)

    analyzer = AnalysisVLM00(5, 10.0)
    analyzer.AddWing(wing)
    analyzer.Analyze([100, 0, 0])

    wingPoints = np.zeros([len(analyzer.wingPanels), 3])
    wingNormals = np.zeros(wingPoints.shape)
    wakePoints = np.zeros([len(analyzer.wakePanels), 3])
    wakeNormals = np.zeros(wakePoints.shape)

    for i in range(0, len(analyzer.wingPanels)):
        wingPoints[i] = analyzer.wingPanels[i][AnalysisIndex.CollacationPoint.value]
        wingNormals[i] = analyzer.wingPanels[i][AnalysisIndex.NormalVector.value]

    for i in range(0, len(analyzer.wakePanels)):
        wakePoints[i] = analyzer.wakePanels[i][AnalysisIndex.CollacationPoint.value]
        wakeNormals[i] = analyzer.wakePanels[i][AnalysisIndex.NormalVector.value]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(wakePoints[:, 0], wakePoints[:, 1], wakePoints[:, 2], c='r')
    ax.scatter(wingPoints[:, 0], wingPoints[:, 1], wingPoints[:, 2], c='g')
    ax.quiver(wakePoints[:, 0], wakePoints[:, 1], wakePoints[:, 2],
              wakeNormals[:, 0], wakeNormals[:, 1], wakeNormals[:, 2], length=0.3, pivot='tail')
    ax.quiver(wingPoints[:, 0], wingPoints[:, 1], wingPoints[:, 2],
              wingNormals[:, 0], wingNormals[:, 1], wingNormals[:, 2], length=0.3, pivot='tail')
    aeroutil.set3DAxesEqual(ax, np.append(wakePoints[:, 0], wingPoints[:, 0]),
                                np.append(wakePoints[:, 1], wingPoints[:, 1]),
                                np.append(wakePoints[:, 2], wingPoints[:, 2]))

def __testVLM00UnitAxes__():
    foil = aerofoil4.AirfoilNACA4Series('3412', aerogen.GeneratorDoubleChebyshev(0, 1, 15))
    wing = aerowing.Wing(foil, aerogen.GeneratorPiecewiseLinear([0.6, 2.0, 0.6], [5, 5]),
        aerogen.GeneratorLinear(0, 10, 10), 0.0, 0, [0, -5, 0])

    analyzer = AnalysisVLM00(8, 15.0)
    analyzer.AddWing(wing)
    analyzer.Analyze([100, 0, 0])

    wingPoints = np.zeros([len(analyzer.wingPanels), 3])
    wingXAxis = np.zeros(wingPoints.shape)
    wingYAxis = np.zeros(wingPoints.shape)
    wingZAxis = np.zeros(wingPoints.shape)

    for i in range(0, len(analyzer.wingPanels)):
        wingPoints[i] = analyzer.wingPanels[i][AnalysisIndex.CollacationPoint.value]
        wingXAxis[i] = analyzer.wingPanels[i][AnalysisIndex.XVector.value]
        wingYAxis[i] = analyzer.wingPanels[i][AnalysisIndex.YVector.value]
        wingZAxis[i] = analyzer.wingPanels[i][AnalysisIndex.NormalVector.value]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(wingPoints[:, 0], wingPoints[:, 1], wingPoints[:, 2], c='k')
    ax.quiver(wingPoints[:, 0], wingPoints[:, 1], wingPoints[:, 2],
              wingXAxis[:, 0], wingXAxis[:, 1], wingXAxis[:, 2], color='r', length=0.1)
    ax.quiver(wingPoints[:, 0], wingPoints[:, 1], wingPoints[:, 2],
              wingYAxis[:, 0], wingYAxis[:, 1], wingYAxis[:, 2], color='g', length=0.1)
    ax.quiver(wingPoints[:, 0], wingPoints[:, 1], wingPoints[:, 2],
              wingZAxis[:, 0], wingZAxis[:, 1], wingZAxis[:, 2], color='b', length=0.1)
    aeroutil.set3DAxesEqual(ax, wingPoints[:, 0], wingPoints[:, 1], wingPoints[:, 2])

def __testVLM00Parameters__():
    foil = aerofoil4.AirfoilNACA4Series('0036', aerogen.GeneratorDoubleChebyshev(0, 1, 15))
    wing = aerowing.Wing(foil, aerogen.GeneratorPiecewiseLinear([0.3, 1.0, 0.3], [5, 6]),
        aerogen.GeneratorLinear(0, 5, 11), 0, 0, [0, -5, 0])

    analyzer = AnalysisVLM00(8, 15.0)
    analyzer.AddWing(wing)
    analyzer.Analyze([50, 0, 0])

    upper, lower = analyzer.GetParameters(0)

    cmap = plt.get_cmap('jet')

    values = np.zeros(2 * len(upper))
    color = np.zeros([2 * len(upper), 4])

    for i in range(0, len(upper)):
        values[2 * i] = upper[i][AnalysisIndex.Velocity.value][2]
        values[2 * i + 1] = values[2 * i]

    minVal = min(values)
    maxVal = max(values)
    for i in range(0, len(values)):
        color[i] = cmap((values[i] - minVal) / (maxVal - minVal))

    print('min value:', minVal)
    print('max value:', maxVal)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.plot_trisurf(wing.upperCoordinates[:, 0], wing.upperCoordinates[:, 1],
    #                wing.upperCoordinates[:, 2], triangles=wing.upperTriangulation,
    #                facecolor=color)
    aeroutil.plot3DColors(ax, wing.upperCoordinates[:, 0], wing.upperCoordinates[:, 1],
                          wing.upperCoordinates[:, 2], wing.upperTriangulation, color)
    aeroutil.set3DAxesEqual(ax, wing.upperCoordinates[:, 0],
                            wing.upperCoordinates[:, 1], wing.upperCoordinates[:, 2])

def __testVLM00Velocity__():
    foil = aerofoil4.AirfoilNACA4Series('3412', aerogen.GeneratorDoubleChebyshev(0, 1, 20))
    wing = aerowing.Wing(foil, aerogen.GeneratorPiecewiseLinear([0.6, 2.0, 0.6], [8, 8]),
        aerogen.GeneratorLinear(0, 10, 16), 0.0, 0, [0, -5, 0])

    analyzer = AnalysisVLM00(8, 15.0)
    analyzer.AddWing(wing)
    analyzer.Analyze([100, 0, 0])

    wingPoints = np.zeros([len(analyzer.wingPanels), 3])
    wingNormals = np.zeros(wingPoints.shape)
    wingVelocity = np.zeros(wingPoints.shape)
    wingVelocityMagnitude = np.zeros(wingPoints.shape[0])
    wakePoints = np.zeros([len(analyzer.wakePanels), 3])
    wakeNormals = np.zeros(wakePoints.shape)

    for i in range(0, len(analyzer.wingPanels)):
        wingPoints[i] = analyzer.wingPanels[i][AnalysisIndex.CollacationPoint.value]
        wingNormals[i] = analyzer.wingPanels[i][AnalysisIndex.NormalVector.value]
        wingVelocity[i] = analyzer.wingPanels[i][AnalysisIndex.Velocity.value]
        wingVelocityMagnitude[i] = np.linalg.norm(wingVelocity[i])

    for i in range(0, len(analyzer.wakePanels)):
        wakePoints[i] = analyzer.wakePanels[i][AnalysisIndex.CollacationPoint.value]
        wakeNormals[i] = analyzer.wakePanels[i][AnalysisIndex.NormalVector.value]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(wakePoints[:, 0], wakePoints[:, 1], wakePoints[:, 2], c='r')
    ax.scatter(wingPoints[:, 0], wingPoints[:, 1], wingPoints[:, 2], c='g')
    ax.quiver(wingPoints[:, 0], wingPoints[:, 1], wingPoints[:, 2],
              wingVelocity[:, 0], wingVelocity[:, 1], wingPoints[:, 2],
              length=0.3, pivot='tail')
    aeroutil.set3DAxesEqual(ax, wingPoints[:, 0], wingPoints[:, 1], wingPoints[:, 2])

#__testVLM00Geometry__()
#__testVLM00UnitAxes__()
__testVLM00Parameters__()
#__testVLM00Velocity__()
