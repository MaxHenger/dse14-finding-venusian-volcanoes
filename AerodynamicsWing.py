# -*- coding: utf-8 -*-
"""
Created on Mon May 16 23:56:09 2016

@author: MaxHenger
"""

import AerodynamicsAirfoil as aerofoil
import AerodynamicsUtil as aeroutil
import AerodynamicsGenerator as aerogen
import AerodynamicsAirfoilNACA4Series as aerofoil4

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Wing:
    def __init__(self, airfoil, chord, halfSpan, sweep, dihedral,
                 translation=[0, 0, 0], rotation=0):
        # Check input for validity
        if not isinstance(airfoil, aerofoil.AirfoilGenerator):
            raise ValueError("Expected 'airfoil' to subclass 'AirfoilGenerator' " +
                "in the 'Airfoil' file")

        if dihedral < -np.pi / 2.0 or dihedral > np.pi / 2.0:
            raise ValueError("Expected 'dihedral' to be within -pi / 2 and pi / 2")

        if isinstance(chord, aerogen.Generator):
            chord = chord.Get()
        elif not aeroutil.isArray(chord):
            raise ValueError("Variable 'chord' should be either an array or " +
                "an instance of AerodynamicsGenerator.Generator")

        if aeroutil.isArray(halfSpan):
            if not aeroutil.isAscending(halfSpan):
                raise ValueError("Expected variable 'halfSpan' to contain only " +
                    "ascending values when it is an array")

            if min(halfSpan) < 0.0:
                raise ValueError("Expected variable 'halfSpan' to contain only " +
                    "numbers larger than 0.0 when it is an array")
        elif isinstance(halfSpan, aerogen.Generator):
            halfSpan = halfSpan.Get()
        else:
            raise ValueError("Variable 'halfSpan' should be either an array or " +
                "an instance of AerodynamicsGenerator.Generator")

        if len(chord) != len(halfSpan):
            raise ValueError("Expected the length of the variable 'chord' to match " +
                "the length of the variable 'halfSpan'")

        if isinstance(sweep, aerogen.Generator):
            sweep = sweep.Get()
        elif isinstance(sweep, float) or isinstance(sweep, int):
            sweep = np.repeat(sweep, len(chord))
        elif not aeroutil.isArray(sweep):
            raise ValueError("Variable 'sweep' should be either an array, a " +
                "float, or an instance of AerodynamicsGenerator.Generator")

        if len(chord) != len(sweep):
            raise ValueError("Expected the length of the variable 'chord' to match " +
                "the length of the variable 'sweep'")

        # Retrieve often-used variables
        upperSurface = airfoil.GetUpperSurface()
        lowerSurface = airfoil.GetLowerSurface()

        # Generate half a wing (in the positive y-plane, starting at x = 0), do
        # not apply dihedral yet, but compensate the half span for the dihedral
        # before rotating
        overCosDihedral = 1.0 / (np.cos(dihedral))
        self.upperCoordinates = np.zeros([upperSurface.shape[0] * len(halfSpan), 3])
        self.lowerCoordinates = np.zeros([lowerSurface.shape[0] * len(halfSpan), 3])
        lastChord = chord[0]
        lastBaseX = 0.0
        lastY = 0.0

        for iSpan in range(0, len(halfSpan)):
            # Generate y-based (base) coordinates
            curY = halfSpan[iSpan]
            curChord = chord[iSpan]
            curSweep = sweep[iSpan]

            baseX = lastBaseX + 0.25 * (lastChord - curChord) + (curY - lastY) * np.sin(curSweep)

            for iChord in range(0, upperSurface.shape[0]):
                localY = curY * overCosDihedral
                localZ = upperSurface[iChord,1] * curChord

                self.upperCoordinates[iSpan * upperSurface.shape[0] + iChord] = [
                    baseX + upperSurface[iChord,0] * curChord,
                    np.cos(dihedral) * localY - np.sin(dihedral) * localZ,
                    np.sin(dihedral) * localY + np.cos(dihedral) * localZ
                ]

            for iChord in range(0, lowerSurface.shape[0]):
                localY = curY * overCosDihedral
                localZ = lowerSurface[iChord,1] * curChord

                self.lowerCoordinates[iSpan * lowerSurface.shape[0] + iChord] = [
                    baseX + lowerSurface[iChord,0] * curChord,
                    np.cos(dihedral) * localY - np.sin(dihedral) * localZ,
                    np.sin(dihedral) * localY + np.cos(dihedral) * localZ
                ]

            lastChord = curChord
            lastBaseX = baseX
            lastY = curY

        # For the first set of coordinates interpolate the coordinates to the
        # initial wingspan point such that it forms a neat axis-aligned plane
        for iChord in range(0, upperSurface.shape[0]):
            point1 = self.upperCoordinates[iChord + upperSurface.shape[0]]
            direction = self.upperCoordinates[iChord] - point1
            factor = (halfSpan[0] - point1[1]) / direction[1]

            if factor < 0:
                raise RuntimeError("Wingspan is ordered such that aligning the " +
                    "first upper chordwise cross-section with the xz-plane results " +
                    "in an overlap. This algorithm cannot deal with this (yet)")

            self.upperCoordinates[iChord] = point1 + factor * direction

        for iChord in range(0, lowerSurface.shape[0]):
            point1 = self.lowerCoordinates[iChord + lowerSurface.shape[0]]
            direction = self.lowerCoordinates[iChord] - point1
            factor = (halfSpan[0] - point1[1]) / direction[1]

            if factor < 0:
                raise RuntimeError("Wingspan is ordered such that aligning the " +
                    "first lower chordwise cross-section with the xz-plane results " +
                    "in an overlap. This algorithm cannot deal with this (yet)")
            self.lowerCoordinates[iChord] = point1 + factor * direction

        # Apply post-generation rotation around the x-axis and consecutive
        # translation
        for i in range(0, len(self.upperCoordinates)):
            curY = self.upperCoordinates[i, 1]
            curZ = self.upperCoordinates[i, 2]
            self.upperCoordinates[i, 0] += translation[0]
            self.upperCoordinates[i, 1] = np.cos(rotation) * curY - np.sin(rotation) * curZ + translation[1]
            self.upperCoordinates[i, 2] = np.sin(rotation) * curY + np.cos(rotation) * curZ + translation[2]

        for i in range(0, len(self.lowerCoordinates)):
            curY = self.lowerCoordinates[i, 1]
            curZ = self.lowerCoordinates[i, 2]
            self.lowerCoordinates[i, 0] += translation[0]
            self.lowerCoordinates[i, 1] = np.cos(rotation) * curY - np.sin(rotation) * curZ + translation[1]
            self.lowerCoordinates[i, 2] = np.sin(rotation) * curY + np.cos(rotation) * curZ + translation[2]

        # Construct triangulations for easy plotting.
        numQuadsX = upperSurface.shape[0] - 1
        numQuadsY = len(halfSpan) - 1
        numUpperQuads = numQuadsX * numQuadsY
        self.upperTriangulation = np.zeros([2 * numUpperQuads, 3])
        self.upperIndices = np.zeros([numQuadsX, numQuadsY, 4], dtype=np.int32)
        self.upperArea = np.zeros([numQuadsX, numQuadsY])
        self.upperCollacation = np.zeros([numQuadsX, numQuadsY, 3])
        self.upperNormal = np.zeros([numQuadsX, numQuadsY, 3])

        for iSpan in range(0, numQuadsY):
            for iChord in range(0, numQuadsX):
                # Part calculating the triangulation indices
                iTriangleBase = 2 * (iSpan * (upperSurface.shape[0] - 1) + iChord)
                iCoordinateBase = iSpan * upperSurface.shape[0] + iChord

                # Determine the corner indices and assemble the triangulation
                br = iCoordinateBase
                tr = iCoordinateBase + 1
                tl = iCoordinateBase + upperSurface.shape[0] + 1
                bl = iCoordinateBase + upperSurface.shape[0]
                self.upperIndices[iChord, iSpan] = [ br, tr, tl, bl ]

                self.upperTriangulation[iTriangleBase] = [ br, tr, bl ]
                self.upperTriangulation[iTriangleBase + 1] = [ tr, tl, bl ]

                # Calculate the collacation points, the surface normals and the
                # surface area
                self.upperCollacation[iChord, iSpan] = (
                    self.upperCoordinates[br] + self.upperCoordinates[tr] +
                    self.upperCoordinates[tl] + self.upperCoordinates[bl]
                ) / 4.0

                toTL = self.upperCoordinates[tl] - self.upperCoordinates[bl]
                toTR = self.upperCoordinates[tr] - self.upperCoordinates[bl]
                toBR = self.upperCoordinates[br] - self.upperCoordinates[bl]

                crossTLAndTR = np.cross(toTL, toTR)
                crossTRAndBR = np.cross(toTR, toBR)

                self.upperArea[iChord, iSpan] = 0.5 * \
                    (np.linalg.norm(crossTLAndTR) + np.linalg.norm(crossTRAndBR))

                crossCombined = crossTLAndTR + crossTRAndBR
                self.upperNormal[iChord, iSpan] = -crossCombined / np.linalg.norm(crossCombined)

        numQuadsX = lowerSurface.shape[0] - 1
        numLowerQuads = numQuadsX * numQuadsY
        self.lowerTriangulation = np.zeros([2 * numLowerQuads, 3])
        self.lowerIndices = np.zeros([numQuadsX, numQuadsY, 4], dtype=np.int32)
        self.lowerArea = np.zeros([numQuadsX, numQuadsY])
        self.lowerCollacation = np.zeros([numQuadsX, numQuadsY, 3])
        self.lowerNormal = np.zeros([numQuadsX, numQuadsY, 3])

        for iSpan in range(0, len(halfSpan) - 1):
            for iChord in range(0, lowerSurface.shape[0] - 1):
                # Part calculating the triangulation indices
                iTriangleBase = 2 * (iSpan * (lowerSurface.shape[0] -1) + iChord)
                iCoordinateBase = iSpan * lowerSurface.shape[0] + iChord

                # Determine the corder indices
                br = iCoordinateBase
                tr = iCoordinateBase + 1
                tl = iCoordinateBase + lowerSurface.shape[0] + 1
                bl = iCoordinateBase + lowerSurface.shape[0]
                self.lowerIndices[iChord, iSpan] = [ tr, br, bl, tl ]

                self.lowerTriangulation[iTriangleBase] = [ br, bl, tr ]
                self.lowerTriangulation[iTriangleBase + 1] = [ bl, tl, tr ]

                # Calculate the collocation points, the surface normals and the
                # surface area
                self.lowerCollacation[iChord, iSpan] = (
                    self.lowerCoordinates[br] + self.lowerCoordinates[tr] +
                    self.lowerCoordinates[tl] + self.lowerCoordinates[bl]
                ) / 4.0

                toTL = self.lowerCoordinates[tl] - self.lowerCoordinates[bl]
                toTR = self.lowerCoordinates[tr] - self.lowerCoordinates[bl]
                toBR = self.lowerCoordinates[br] - self.lowerCoordinates[bl]

                crossTLAndTR = np.cross(toTL, toTR)
                crossTRAndBR = np.cross(toTR, toBR)

                self.lowerArea[iChord, iSpan] = 0.5 * \
                    (np.linalg.norm(crossTLAndTR) + np.linalg.norm(crossTRAndBR))

                crossCombined = crossTLAndTR + crossTRAndBR
                self.lowerNormal[iChord, iSpan] = crossCombined / np.linalg.norm(crossCombined)

    # Functions related to lower panels
    def GetNumLowerPanelsX(self):
        return self.lowerArea.shape[0]

    def GetNumLowerPanelsY(self):
        return self.lowerArea.shape[1]

    def GetLowerPanelArea(self, iX, iY):
        return self.lowerArea[iX, iY]

    def GetLowerPanelPoints(self, iX, iY):
        return np.asarray([self.lowerCoordinates[self.lowerIndices[iX, iY, 0]],
                           self.lowerCoordinates[self.lowerIndices[iX, iY, 1]],
                           self.lowerCoordinates[self.lowerIndices[iX, iY, 2]],
                           self.lowerCoordinates[self.lowerIndices[iX, iY, 3]]])

    def GetLowerPanelCollacationPoint(self, iX, iY):
        return self.lowerCollacation[iX, iY]

    def GetLowerPanelNormal(self, iX, iY):
        return self.lowerNormal[iX, iY]

    # Functions related to upper panels
    def GetNumUpperPanelsX(self):
        return self.upperArea.shape[0]

    def GetNumUpperPanelsY(self):
        return self.upperArea.shape[1]

    def GetUpperPanelArea(self, iX, iY):
        return self.upperArea[iX, iY]

    def GetUpperPanelPoints(self, iX, iY):
        return np.asarray([self.upperCoordinates[self.upperIndices[iX, iY, 0]],
                           self.upperCoordinates[self.upperIndices[iX, iY, 1]],
                           self.upperCoordinates[self.upperIndices[iX, iY, 2]],
                           self.upperCoordinates[self.upperIndices[iX, iY, 3]]])

    def GetUpperPanelCollacationPoint(self, iX, iY):
        return self.upperCollacation[iX, iY]

    def GetUpperPanelNormal(self, iX, iY):
        return self.upperNormal[iX, iY]

def __testWing__():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    airfoil = aerofoil4.AirfoilNACA4Series('2436', aerogen.GeneratorDoubleChebyshev(0, 1, 15))
    w = Wing(airfoil, aerogen.GeneratorPiecewiseLinear([3.0, 1.0, 0.2], [3, 7]),
             aerogen.GeneratorLinear(0, 5, 10),
             np.pi / 4.0, np.pi / 6.0)

    ax.plot_trisurf(w.upperCoordinates[:, 0],
                    w.upperCoordinates[:, 1],
                    w.upperCoordinates[:, 2],
                    triangles=w.upperTriangulation)
    ax.plot_trisurf(w.lowerCoordinates[:, 0],
                    w.lowerCoordinates[:, 1],
                    w.lowerCoordinates[:, 2],
                    triangles=w.lowerTriangulation)

    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim( 0.0, 5.0)
    ax.set_zlim(-2.5, 2.5)

def __testWingRotatedTranslated__():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    airfoil = aerofoil4.AirfoilNACA4Series('0030', aerogen.GeneratorDoubleChebyshev(0, 1, 25))
    w = Wing(airfoil, aerogen.GeneratorEllipticLinear(2.0, 0.2, 10),
             aerogen.GeneratorLinear(0, 5, 10), np.pi / 25, 0, [-2, 0, 0], np.pi / 2.0)

    ax.plot_trisurf(w.upperCoordinates[:, 0],
                    w.upperCoordinates[:, 1],
                    w.upperCoordinates[:, 2],
                    triangles=w.upperTriangulation)
    ax.plot_trisurf(w.lowerCoordinates[:, 0],
                    w.lowerCoordinates[:, 1],
                    w.lowerCoordinates[:, 2],
                    triangles=w.lowerTriangulation)

    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.5)
    ax.set_zlim( 0.0, 5.0)

def __testWingNormals__():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    airfoil = aerofoil4.AirfoilNACA4Series('8426', aerogen.GeneratorDoubleChebyshev(0, 1, 15))
    w = Wing(airfoil, aerogen.GeneratorLinear(2.0, 0.6, 10),
             aerogen.GeneratorLinear(0, 4, 10), np.pi / 12.0, 0,
             [-1.0, -1.0, 1.0], -np.pi / 4)

    ax.plot_trisurf(w.upperCoordinates[:, 0],
                    w.upperCoordinates[:, 1],
                    w.upperCoordinates[:, 2],
                    triangles=w.upperTriangulation)
    ax.plot_trisurf(w.lowerCoordinates[:, 0],
                    w.lowerCoordinates[:, 1],
                    w.lowerCoordinates[:, 2],
                    triangles=w.lowerTriangulation)

    # assemble points and vectors of collacation points
    numX = w.GetNumLowerPanelsX()
    numY = w.GetNumLowerPanelsY()

    lowerP = np.zeros([numX * numY, 3])
    lowerN = np.zeros(lowerP.shape)

    for iX in range(0, numX):
        for iY in range(0, numY):
            lowerP[iX * numY + iY] = w.GetLowerPanelCollacationPoint(iX, iY)
            lowerN[iX * numY + iY] = w.GetLowerPanelNormal(iX, iY)

    numX = w.GetNumUpperPanelsX()
    numY = w.GetNumUpperPanelsY()

    upperP = np.zeros([numX * numY, 3])
    upperN = np.zeros(upperP.shape)

    for iX in range(0, numX):
        for iY in range(0, numY):
            upperP[iX * numY + iY] = w.GetUpperPanelCollacationPoint(iX, iY)
            upperN[iX * numY + iY] = w.GetUpperPanelNormal(iX, iY)

    # plot normals using a quiver
    ax.quiver(upperP[:, 0], upperP[:, 1], upperP[:, 2],
              upperN[:, 0], upperN[:, 1], upperN[:, 2],
              length=0.5, pivot='tail', color='g')
    ax.quiver(lowerP[:, 0], lowerP[:, 1], lowerP[:, 2],
              lowerN[:, 0], lowerN[:, 1], lowerN[:, 2],
              length=0.5, pivot='tail', color='r')

    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim( 0.0, 5.0)
    ax.set_zlim(-2.5, 2.5)

def __testWingCenter__():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    airfoil = aerofoil4.AirfoilNACA4Series('8440', aerogen.GeneratorDoubleChebyshev(0, 1, 15))
    w = Wing(airfoil, aerogen.GeneratorLinear(2.0, 0.6, 10),
             aerogen.GeneratorLinear(0, 4, 10), np.pi / 12.0, 0,
             [-1.0, -1.0, 1.0], -np.pi / 4)

    ax.plot_trisurf(w.upperCoordinates[:, 0],
                    w.upperCoordinates[:, 1],
                    w.upperCoordinates[:, 2],
                    triangles=w.upperTriangulation)
    ax.plot_trisurf(w.lowerCoordinates[:, 0],
                    w.lowerCoordinates[:, 1],
                    w.lowerCoordinates[:, 2],
                    triangles=w.lowerTriangulation)

    # Assemble collacation points in a single array
    numX = w.GetNumLowerPanelsX()
    numY = w.GetNumLowerPanelsY()

    upperP = np.zeros([numX * numY, 3])

    for iX in range(0, numX):
        for iY in range(0, numY):
            upperP[iX * numY + iY] = w.GetLowerPanelCollacationPoint(iX, iY)

    numX = w.GetNumUpperPanelsX()
    numY = w.GetNumUpperPanelsY()

    lowerP = np.zeros([numX * numY, 3])

    for iX in range(0, numX):
        for iY in range(0, numY):
            lowerP[iX * numY + iY] = w.GetUpperPanelCollacationPoint(iX, iY)

    ax.scatter(lowerP[:, 0], lowerP[:, 1], lowerP[:, 2], c='r')
    ax.scatter(upperP[:, 0], upperP[:, 1], upperP[:, 2], c='g')

    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim( 0.0, 5.0)
    ax.set_zlim(-2.5, 2.5)

#__testWing__()
#__testWingRotatedTranslated__()
#__testWingNormals__()
#__testWingCenter__()
