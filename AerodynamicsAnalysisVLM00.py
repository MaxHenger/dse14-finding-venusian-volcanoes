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
    NormalVector = 2
    SourceStrength = 3
    DoubletStrength = 4
    InfluenceWingDoublet = 5
    InfluenceWingSource = 6
    InfluenceWakeDoublet = 7
    LocalU = 8
    LocalV = 9
    LocalW = 10
    IsUpper = 11
    WakePanel = 12
    Total = 13

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
        self.Vinf = Vinf

        if self.wingsChanged == True:
            self.__updateWingGeometry__()
            self.__updateInfluenceParameters__()

        self.__updateSingularityStrengths__()

    def __updateWingGeometry__(self):
        self.wingPanels = []
        self.wakePanels = []

        for iWing in range(0, len(self.wings)):
            localWingPanels, localWakePanels = self.__generateWingGeometry__(self.wings[iWing], len(self.wakePanels))
            self.wingPanels.extend(localWingPanels)
            self.wakePanels.extend(localWakePanels)

    def __generateWingGeometry__(self, wing, wakeOffset):
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
                newPanel[AnalysisIndex.CollacationPoint.value] = wing.GetLowerPanelCollacationPoint(iX, iY)
                newPanel[AnalysisIndex.Corners.value] = wing.GetLowerPanelPoints(iX, iY)
                newPanel[AnalysisIndex.NormalVector.value] = wing.GetLowerPanelNormal(iX, iY)
                newPanel[AnalysisIndex.IsUpper.value] = False
                newPanel[AnalysisIndex.WakePanel.value] = -1
                wingPanels.append(newPanel)

            # Add the last panel, which is connected to the wake
            newPanel = [None] * AnalysisIndex.Total.value
            newPanel[AnalysisIndex.CollacationPoint.value] = wing.GetLowerPanelCollacationPoint(numX - 1, iY)
            newPanel[AnalysisIndex.Corners.value] = wing.GetLowerPanelPoints(numX - 1, iY)
            newPanel[AnalysisIndex.NormalVector.value] = wing.GetLowerPanelNormals(numX - 1, iY)
            newPanel[AnalysisIndex.IsUpper.value] = False
            newPanel[AnalysisIndex.WakePanel.value] = wakeOffset + len(wakePanels)
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

                toTL = newPanel[AnalysisIndex.Corners.value][2] - newPanel[AnalysisIndex.Corners.value][3]
                toTR = newPanel[AnalysisIndex.Corners.value][1] - newPanel[AnalysisIndex.Corners.value][3]
                toBR = newPanel[AnalysisIndex.Corners.value][0] - newPanel[AnalysisIndex.Corners.value][3]

                normal = np.cross(toTR, toTL) + np.cross(toBR, toTR)
                normal /= np.linalg.norm(normal)

                newPanel[AnalysisIndex.NormalVector.value] = normal
                wakePanels.append(newPanel)

        numX = wing.GetNumUpperPanelsX()
        numY = wing.GetNumUpperPanelsY()

        for iY in range(0, numY):
            for iX in range(0, numX - 1):
                # Generate new upper wing panel
                newPanel = [None] * AnalysisIndex.Total.value
                newPanel[AnalysisIndex.CollacationPoint.value] = wing.GetUpperPanelCollacationPoint(iX, iY)
                newPanel[AnalysisIndex.Corners.value] = wing.GetUpperPanelPoints(iX, iY)
                newPanel[AnalysisIndex.NormalVector.value] = wing.GetUpperPanelNormal(iX, iY)
                newPanel[AnalysisIndex.IsUpper.value] = True
                newPanel[AnalysisIndex.WakePanel.value] = -1
                wingPanels.append(newPanel)

            # Generate last upper panel connected to the wake
            newPanel = [None] * AnalysisIndex.Total.value
            newPanel[AnalysisIndex.CollacationPoint.value] = wing.GetUpperPanelCollacationPoint(numX - 1, iY)
            newPanel[AnalysisIndex.Corners.value] = wing.GetUpperPanelPoints(numX - 1, iY)
            newPanel[AnalysisIndex.NormalVector.value] = wing.GetUpperPanelNormal(numX - 1, iY)
            newPanel[AnalysisIndex.IsUpper.value] = True
            newPanel[AnalysisIndex.WakePanel.value] = baseWake[iY]
            wingPanels.append(newPanel)

        return wingPanels, wakePanels

    def __updateInfluenceParameters__(self):
        # For each wing panel calculate the influences of all other panels
        for iWingConsidered in range(0, len(self.wingPanels)):
            # Determine influence of other wing panel source and doublet
            # singularity functions
            curPanel = self.wingPanels[iWingConsidered]
            curPanel[AnalysisIndex.InfluenceWingSource.value] = np.zeros(len(self.wingPanels))
            curPanel[AnalysisIndex.InfluenceWingDoublet.value] = np.zeros(len(self.wingPanels))

            for iWingAffecting in range(0, len(self.wingPanels)):
                curPanel[AnalysisIndex.InfluenceWingSource.value][iWingAffecting] = \
                    self.__calculateSourcePanelInfluence__(
                        curPanel[AnalysisIndex.CollacationPoint.value],
                        curPanel[AnalysisIndex.Corners.value],
                        curPanel[AnalysisIndex.NormalVector.value],
                        self.wingPanels[iWingAffecting][AnalysisIndex.CollacationPoint.value]
                    )

                curPanel[AnalysisIndex.InfluenceWingDoublet.value][iWingAffecting] = \
                    self.__calculateDoubletPanelInfluence__(
                        curPanel[AnalysisIndex.CollacationPoint.value],
                        curPanel[AnalysisIndex.Corners.value],
                        curPanel[AnalysisIndex.NormalVector.value],
                        self.wingPanels[iWingAffecting][AnalysisIndex.CollacationPoint.value]
                    )

            # Determine influence of wake panels
            curPanel[AnalysisIndex.InfluenceWakeDoublet.value] = np.zeros(len(self.wakePanels))

            for iWakeAffecting in range(0, len(self.wakePanels)):
                curPanel[AnalysisIndex.InfluenceWakeDoublet.value][iWakeAffecting] = \
                    self.__calculateDoubletPanelInfluence__(
                        curPanel[AnalysisIndex.CollacationPoint.value],
                        curPanel[AnalysisIndex.Corners.value],
                        curPanel[AnalysisIndex.NormalVector.value],
                        self.wakePanels[iWakeAffecting][AnalysisIndex.CollacationPoint.value]
                    )

    def __updateSingularityStrengths__(self, Vinf):
        # Update all wing panel source strengths
        for iWing in range(0, len(self.wingPanels)):
            self.wingPanels[iWing][AnalysisIndex.SourceStrength.value] = \
                np.dot(self.wingPanels[iWing][AnalysisIndex.NormalVector.value], Vinf)

        # Generate the matrix equation that will be solved
        matrix = np.zeros([len(self.wingPanels), len(self.wingPanels)])
        vector = np.zeros(len(self.wingPanels))
        wakeLookup = {}

        for iWingConsidered in range(0, len(self.wingPanels)):
            curWing = self.wingPanels[iWingConsidered]

            for iWingAffecting in range(0, len(self.WingPanels)):
                curAffecting = self.wingPanels[iWingAffecting]

                vector[iWingConsidered] -= \
                    curWing[AnalysisIndex.InfluenceWingSource.value][iWingAffecting] * \
                    curAffecting[AnalysisIndex.SourceStrength.value]

                if curAffecting[AnalysisIndex.WakePanel.value] != -1:
                    # A wake is connected to this panel
                    if self.wingPanels[iWing][AnalysisIndex.IsUpper] == True:
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

    def __calculateAuxilliary__(self, panelCenter, panelCoordinates, panelNormal, point):
        # Note that this algorithm is taken from "Low-Speed Aerodynamics: From
        # Wing Theory to Panel Methods" by Joseph Katz and Allen Plotkin

        # Move all coordinates to the local reference frame
        panelCoordinates[0] -= panelCenter
        panelCoordinates[1] -= panelCenter
        panelCoordinates[2] -= panelCenter
        panelCoordinates[3] -= panelCenter
        point -= panelCenter

        # Create a local x-vector and a local y-vector
        localX = (panelCoordinates[0] + panelCoordinates[1]) / 2.0
        localY = -np.cross(localX, panelNormal)

        # Project the point onto the axes to get a relative vector with respect
        # to the panel's reference frame
        localPoint = [
            np.linalg.dot(point, localX),
            np.linalg.dot(point, localY),
            np.linalg.dot(point, panelNormal)
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
            e[i] = np.power(localPoint[0] - panelCoordinates[i, 0], 2.0) + np.power(localPoint[2], 2.0)
            h[i] = (localPoint[0] - panelCoordinates[i, 0]) * (localPoint[1] - panelCoordinates[i, 1])

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

    def __calculateSourcePanelInfluence__(self, panelCenter, panelCoordinates, panelNormal, point):
        # Note that this algorithm is taken from "Low-Speed Aerodynamics: From
        # Wing Theory to Panel Methods" by Joseph Katz and Allen Plotkin

        # Retrieve auxilliary variables
        (panelCoordinates, localPoint, r, e, h,
         m12, m23, m34, m41, d12, d23, d34, d41) = \
            self.__calculateAuxilliary__(panelCenter, panelCoordinates, panelNormal, point)

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
            # Local u-component of velocity
            ( 1.0 / (4.0 * np.pi) * (
                (panelCoordinates[1, 1] - panelCoordinates[0, 1]) / d12 *
                np.log( (r[0] + r[1] - d12) / (r[0] + r[1] + d12) ) +
                (panelCoordinates[2, 1] - panelCoordinates[1, 1]) / d23 *
                np.log( (r[1] + r[2] - d23) / (r[1] + r[2] + d23) ) +
                (panelCoordinates[3, 1] - panelCoordinates[2, 1]) / d34 *
                np.log( (r[2] + r[3] - d34) / (r[2] + r[3] + d34) ) +
                (panelCoordinates[0, 1] - panelCoordinates[3, 1]) / d41 *
                np.log( (r[3] + r[0] - d41) / (r[3] + r[0] + d41) )
            )),
            # Local v-component of velocity
            ( 1.0 / (4.0 * np.pi) * (
                (panelCoordinates[0, 0] - panelCoordinates[1, 0]) / d12 *
                np.log( (r[0] + r[1] - d12) / (r[0] + r[1] + d12) ) +
                (panelCoordinates[1, 0] - panelCoordinates[2, 0]) / d23 *
                np.log( (r[1] + r[2] - d23) / (r[1] + r[2] + d23) ) +
                (panelCoordinates[2, 0] - panelCoordinates[3, 0]) / d34 *
                np.log( (r[2] + r[3] - d34) / (r[2] + r[3] + d34) ) +
                (panelCoordinates[3, 0] - panelCoordinates[0, 0]) / d41 *
                np.log( (r[3] + r[0] - d41) / (r[3] + r[0] + d41) )
            )),
            # Local w-component of velocity
            ( 1.0 / (4.0 * np.pi) * (
                np.arctan2(m12 * e[0] - h[0], localPoint[2] * r[0]) -
                np.arctan2(m12 * e[1] - h[1], localPoint[2] * r[1]) +
                np.arctan2(m23 * e[1] - h[1], localPoint[2] * r[1]) -
                np.arctan2(m23 * e[2] - h[2], localPoint[2] * r[2]) +
                np.arctan2(m34 * e[2] - h[2], localPoint[2] * r[2]) -
                np.arctan2(m34 * e[3] - h[3], localPoint[2] * r[3]) +
                np.arctan2(m41 * e[3] - h[3], localPoint[2] * r[3]) -
                np.arctan2(m41 * e[0] - h[0], localPoint[2] * r[0])
            ))
        )

    def __calculateDoubletPanelInfluence__(self, panelCenter, panelCoordinates, panelNormal, point):
        # Note that this algorithm is taken from "Low-Speed Aerodynamics: From
        # Wing Theory to Panel Methods" by Joseph Katz and Allen Plotkin

        # Retrieve auxilliary variables
        (_, localPoint, r, e, h, m12, m23, m34, m41, d12, d23, d34, d41) = \
            self.__calculateAuxilliary__(panelCenter, panelCoordinates, panelNormal, point)

        # Calculate total influence
        return \
        (   # Influence parameters
            1.0 / (4.0 * np.pi) * (
                np.arctan2(m12 * e[0] * h[0], localPoint[2] * r[0]) -
                np.arctan2(m12 * e[1] - h[1], localPoint[2] * r[1]) +
                np.arctan2(m23 * e[1] - h[1], localPoint[2] * r[1]) -
                np.arctan2(m23 * e[2] - h[2], localPoint[2] * r[2]) +
                np.arctan2(m34 * e[2] - h[2], localPoint[2] * r[2]) -
                np.arctan2(m34 * e[3] - h[3], localPoint[2] * r[3]) +
                np.arctan2(m41 * e[4] - h[4], localPoint[2] * r[3]) -
                np.arctan2(m41 * e[0] - h[0], localPoint[2] * r[0])
            ),
            # Local u-component of velocity
            1.0 / (4.0 * np.pi) * (
                localPoint[2] * (panelCoordinates[0, 1] - panelCoordinates[1, 1]) * (r[0] + r[1]) /
                (
                    r[0] * r[1] * (r[1] * r[2] - (
                        (localPoint[0] - panelCoordinates[0, 0]) *
                        (localPoint[0] - panelCoordinates[1, 0]) +
                        (localPoint[1] - panelCoordinates[0, 1]) *
                        (localPoint[1] - panelCoordinates[1, 1]) +
                        np.power(localPoint[2], 2.0)
                    ))
                ) +
                localPoint[2] * (panelCoordinates[1, 1] - panelCoordinates[2, 1]) * (r[1] + r[2]) /
                (
                    # CONTINUE HERE, SEE P. 153 OF DOCUMENT
                )
            )
        )

def __testVLM00Geometry__():
    foil = aerofoil4.AirfoilNACA4Series('2436',
        aerogen.GeneratorSingleChebyshev(0.0, 1.0, 12))

    wing = aerowing.Wing(foil, aerogen.GeneratorLinear(2.0, 0.6, 10),
        aerogen.GeneratorLinear(0, 5, 10), np.pi / 12.0, 0)

    analyzer = AnalysisVLM00(5, 10.0)
    analyzer.AddWing(wing)
    analyzer.Analyze()

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

#__testVLM00Geometry__()
