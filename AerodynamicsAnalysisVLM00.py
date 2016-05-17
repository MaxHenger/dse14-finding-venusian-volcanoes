# -*- coding: utf-8 -*-
"""
Created on Tue May 17 17:02:34 2016

@author: MaxHenger
"""

import AerodynamicsWing as aerowing

import numpy as np
from enum import Enum

class AnalysisIndex(Enum):
    CollacationPoint = 0
    CornerBR = 1
    CornerTR = 2
    CornerTL = 3
    CornerBL = 4
    NormalVector = 5
    SourceStrength = 6
    DoubletStrength = 7
    InfluenceWingDoublet = 8
    InfluenceWingSource = 9
    InfluenceWakeDoublet = 10
    Total = 11

class AnalysisVLM00:
    def __init__(self, numWake, wakeEnd):
        self.wings = []
        self.wingsChanged = True
        self.numWake = numWake
        self.wakeEnd = wakeEnd

    def AddWing(self, wing):
        # Check if the wing is the expected AerodynamicsWing.Wing object
        if not isinstance(wing, aerowing.Wing):
            raise ValueError("Expected variable 'wing' to be an instance of the " +
                "AerodynamicsWing.Wing class")

        self.wings.append(wing)

    def Analyze(self):
        if (self.wingsChanged):
            return True

    def __generateWingGeometry__(self, wing):
        # Convert wing data into a set of upper surface-, lower surface- and
        # wake panels
        numLowerX = wing.GetNumLowerPanelsX()
        numLowerY = wing.GetNumLowerPanelsY()
        numUpperX = wing.GetNumUpperPanelsX()
        numUpperY = wing.GetNumUpperPanelsY()

        wakePanels = []
        wingPanels = []

        numX = w.GetNumLowerPanelsX()
        numY = w.GetNumLowerPanelsY()

        # Move along span, do each row of upper and lower panels while
        # generating the wake. It is assumed that the airfoil closes neatly
        # such that a wake only has to be attached to one wake panel
        for iY in range(0, numY):
            for iX in range(0, numX):
                # Generate new lower wing panel
                newPanel = [None] * AnalysisIndex.Total
                newPanel[AnalysisIndex.CollacationPoint] = wing.GetLowerPanelCollacationPoint(iX, iY)
                newPanel[AnalysisIndex.Corners] = wing.GetLowerPanelPoints(iX, iY)
                newPanel[AnalysisIndex.NormalVector] = wing.GetLowerPanelNormal(iX, iY)
                wingPanels.append(newPanel)

            # Generate new wake panel
            tempCoords = wing.GetLowerPanelPoints(iX, iY)
            pointInner = tempCoords[0]
            pointOuter = tempCoords[1]

            for iW in range(0, self.numWake)


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
        return -1.0 / (4.0 * np.pi) * (
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
        )

    def __calculateDoubletPanelInfluence__(self, panelCenter, panelCoordinates, panelNormal, point):
        # Note that this algorithm is taken from "Low-Speed Aerodynamics: From
        # Wing Theory to Panel Methods" by Joseph Katz and Allen Plotkin

        # Retrieve auxilliary variables
        (_, localPoint, r, e, h, m12, m23, m34, m41, d12, d23, d34, d41) = \
            self.__calculateAuxilliary__(panelCenter, panelCoordinates, panelNormal, point)

        # Calculate total influence
        return 1.0 / (4.0 * np.pi) * (
            np.arctan2(m12 * e[0] * h[0], localPoint[2] * r[0]) -
            np.arctan2(m12 * e[1] - h[1], localPoint[2] * r[1]) +
            np.arctan2(m23 * e[1] - h[1], localPoint[2] * r[1]) -
            np.arctan2(m23 * e[2] - h[2], localPoint[2] * r[2]) +
            np.arctan2(m34 * e[2] - h[2], localPoint[2] * r[2]) -
            np.arctan2(m34 * e[3] - h[3], localPoint[2] * r[3]) +
            np.arctan2(m41 * e[4] - h[4], localPoint[2] * r[3]) -
            np.arctan2(m41 * e[0] - h[0], localPoint[2] * r[0])
        )
