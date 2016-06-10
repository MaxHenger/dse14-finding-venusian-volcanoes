# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 17:03:18 2016

TrackContour defines a class that can be used that generated line-based contours
for a dataset and easily reducing an area to the required data. Note that it
works perfectly as long as the maps are fine enough. I suspect that errors will
occur if the value goes from below vMin to vMax between two datapoints.

According to the old dutch proverb: A warned man counts for two.

@author: MaxHenger
"""

import TrackCommon
import TrackLookup

import numpy as np

class ContourBaseSection:
    def __init__(self, parent, includer=None, includeBlobIndex=None,
                 excluders=None, excludeBlobIndices=None):
        if excluders != None or excludeBlobIndices != None:
            if len(excluders) != len(excludeBlobIndices):
                raise ValueError("Expected 'excluders' length to equal 'includeBlobIndices'")

            self.excluders = excluders
            self.excludeBlobIndices = excludeBlobIndices
        else:
            self.excluders = []
            self.excludeBlobIndices = []

        self.includer = includer
        self.includeBlobIndex = includeBlobIndex
        self.parent = parent

    def isInside(self, x, y):
        # This is the base class
        raise ValueError('isInside called on the ContourBaseSection class. ' +
                         'This is illegal as it is a base class')
        return False

    def setIncluder(self, includer, blobIndex):
        self.includer = includer
        self.includeBlobIndex = blobIndex

    def addExcluder(self, excluder, blobIndex):
        self.excluders.append(excluder)
        self.excludeBlobIndices.append(blobIndex)

    def getIncluderBlobIndex(self):
        return self.includeBlobIndex

    def getIncluder(self):
        return self.includer

    def getNumExcluders(self):
        return len(self.excluders)

    def getExcluderBlobIndex(self, index):
        return self.excludeBlobIndices[index]

    def getExcluder(self, index):
        return self.excluders[index]

class ContourSingleSection(ContourBaseSection):
    def __init__(self, parent, includer=None, includeBlobIndex=None,
                 excluders=None, excludeBlobIndices=None):
        super(ContourSingleSection, self).__init__(parent, includer,
            includeBlobIndex, excluders, excludeBlobIndices)

    def isInside(self, x, y):
        if x < self.parent.axisX[0] or x > self.parent.axisX[-1]:
            raise ValueError('x is outside of dataset bounds')

        if y < self.parent.axisY[0] or y > self.parent.axisY[-1]:
            raise ValueError('y is outside of dataset bounds')

        iX = self.parent.axisXLookup(self.parent.axisX, x)
        iY = self.parent.axisYLookup(self.parent.axisY, y)

        val = 0.0

        if iX == len(self.parent.axisX) - 1:
            if iY == len(self.parent.axisY) - 1:
                if self.parent.blobMap[iY, iX] != self.includeBlobIndex:
                    return False

                if self.parent.data[iY, iX] < self.parent.vMin:
                    return False

                if self.parent.data[iY, iX] > self.parent.vMax:
                    return False

                return True

            # iY is not at the edge
            if self.parent.blobMap[iY, iX] != self.includeBlobIndex and \
                    self.parent.blobMap[iY + 1, iX] != self.includeBlobIndex:
                return False

            val = TrackCommon.Lerp(self.parent.axisY[iY], self.parent.data[iY, iX],
                                   self.parent.axisY[iY + 1], self.parent.data[iY + 1, iX], y)

            if val < self.parent.vMin:
                return False

            if val > self.parent.vMax:
                return False

            return True

        if iY == len(self.parent.axisY) - 1:
            # iX is not at the edge
            if self.parent.blobMap[iY, iX] != self.includeBlobIndex and \
                    self.parent.blobMap[iY, iX + 1] != self.includeBlobIndex:
                return False

            val = TrackCommon.Lerp(self.parent.axisX[iX], self.parent.data[iY, iX],
                                   self.parent.axisX[iX + 1], self.parent.data[iY, iX + 1], x)

            if val < self.parent.vMin:
                return False

            if val > self.parent.vMax:
                return False

            return True

        # Both are not at the edge: perform bilinear interpolation
        if self.parent.blobMap[iY, iX] != self.includeBlobIndex and \
                self.parent.blobMap[iY + 1, iX] != self.includeBlobIndex and \
                self.parent.blobMap[iY, iX + 1] != self.includeBlobIndex and \
                self.parent.blobMap[iY + 1, iX + 1] != self.includeBlobIndex:
            return False

        val = TrackCommon.Bilerp(self.parent.axisX[iX], self.parent.axisX[iX + 1],
                                 self.parent.axisY[iY], self.parent.axisY[iY + 1],
                                 self.parent.data[iY, iX], self.parent.data[iY + 1, iX],
                                 self.parent.data[iY, iX + 1], self.parent.data[iY + 1, iX + 1],
                                 x, y)

        if val < self.parent.vMin:
            return False

        if val > self.parent.vMax:
            return False

        return True

    def _rangeFromIndices_(self, index0, index1, axis, values):
        indices = [index0, index1]
        result = [0, 0]

        for i in range(0, 2):
            vTarget = 0
            iOther = indices[i] + 2 * i - 1

            if indices[i] == 0 or indices[i] == len(axis) - 1:
                result[i] = axis[indices[i]]
            else:
                if values[indices[i]] > self.parent.vMax:
                    vTarget = self.parent.vMax
                else:
                    vTarget = self.parent.vMin

                result[i] = TrackCommon.Lerp(values[indices[i]],
                    axis[indices[i]], values[iOther], axis[iOther], vTarget)

        return result

    def _ranges_(self, axis, values, isBlob):
        ranges = []
        isPart = isBlob[0] and values[0] >= self.parent.vMin and \
            values[0] <= self.parent.vMax
        oldIndex = 0

        for i in range(1, len(values)):
            if isBlob[i] and values[i] >= self.parent.vMin and \
                    values[i] <= self.parent.vMax:
                # Current index points to a part that is inside the valid range
                if not isPart:
                    isPart = True
                    oldIndex = i
            elif isPart:
                # Current index points to a part that is not inside the valid
                # range.
                isPart = False
                ranges.append(self._rangeFromIndices_(oldIndex, i - 1, axis, values))

        if isPart:
            ranges.append(self._rangeFromIndices_(oldIndex, len(values) - 1, axis, values))
            
        return np.asarray(ranges)

    def getHorizontalRanges(self, y):
        if y < self.parent.axisY[0] or y > self.parent.axisY[-1]:
            raise ValueError("y is outside of dataset bounds")

        # find the location of the y-coordinate
        iY = self.parent.axisYLookup(self.parent.axisY, y)

        #print(' * at ', round(y, 5), '( index', iY, ')')

        valueArray = None
        isBlob = np.empty([len(self.parent.axisX)], dtype=np.bool)

        if iY == len(self.parent.axisY) - 1:
            # Special case
            valueArray = self.parent.data[iY, :]

            for i in range(0, len(self.parent.axisX)):
                isBlob[i] = self.parent.blobMap[iY, i] == self.includeBlobIndex
        else:
            valueArray = TrackCommon.Lerp(self.parent.axisY[iY], self.parent.data[iY, :],
                                          self.parent.axisY[iY + 1], self.parent.data[iY + 1, :], y)

            for i in range(0, len(self.parent.axisX)):
                isBlob[i] = self.parent.blobMap[iY, i] == self.includeBlobIndex or \
                    self.parent.blobMap[iY + 1, i] == self.includeBlobIndex

        #print(' > values:', valueArray)

        # Convert the generated value array into a set of ranges which lie
        # inside the target region of the contour
        return self._ranges_(self.parent.axisX, valueArray, isBlob)

    def getVerticalRanges(self, x):
        if x < self.parent.axisX[0] or x > self.parent.axisX[-1]:
            raise ValueError("x is outside of dataset bounds")

        # Find the location of the x-coordinate
        iX = self.parent.axisXLookup(self.parent.axisX, x)

        valueArray = None
        isBlob = np.empty([len(self.parent.axisY)], dtype=np.bool)

        if iX == len(self.parent.axisX) - 1:
            # Special case
            valueArray = self.parent.data[:, iX]

            for i in range(0, len(self.parent.axisY)):
                isBlob[i] = self.parent.blobMap[i, iX] == self.includeBlobIndex
        else:
            valueArray = TrackCommon.Lerp(self.parent.axisX[iX], self.parent.data[:, iX],
                                          self.parent.axisX[iX + 1], self.parent.data[:, iX + 1], x)

            for i in range(0, len(self.parent.axisY)):
                isBlob[i] = self.parent.blobMap[i, iX] == self.includeBlobIndex or \
                    self.parent.blobMap[i, iX + 1] == self.includeBlobIndex

        # Move along the value array and get the indices of the ranges of values
        # that lie inside the contour
        return self._ranges_(self.parent.axisY, valueArray, isBlob)

class Contour:
    # Constants for readability (and screwing up efficiency, damn you python!).
    # The index is important and may not change as the pixel-finding algorithm
    # depends on these values
    LEFT = 0
    TOPLEFT = 1
    TOP = 2
    TOPRIGHT = 3
    RIGHT = 4
    BOTTOMRIGHT = 5
    BOTTOM = 6
    BOTTOMLEFT = 7
    NONE = 8
    START = 3 # See the _nextEdge_ loop for the reason behind this start value

    # Bitmaps used to determine which interpolations to perform when moving from
    # one data point to the other. These values are used in the interpolateDir
    # map
    BITMAP_LEFT = 0b0001
    BITMAP_TOP = 0b0010
    BITMAP_RIGHT = 0b0100
    BITMAP_BOTTOM = 0b1000

    # directionWord is a map to make debugging easier with all these direction
    # numbers being thrown around
    directionWord = ['L ', 'TL', 'T ', 'TR', 'R ', 'BR', 'B ', 'BL', 'ERROR']

    # tileOffsets consists of the often-used offsets to perform pixel edge
    # walking. Each entry consists of a [yOffset, xOffset] pair
    tileOffsets = [
        [0, -1], [1, -1], [1, 0], [1, 1],
        [0, 1], [-1, 1], [-1, 0], [-1, -1]
    ]

    # interpolateDir is a lookup map that is used to determine in which
    # direction to interpolate to find the contour. The first index is the
    # original direction, the second index is the new direction. The bitmap
    # describes the BITMAP_LEFT, etc. values in the first four bits. The second
    # two bits describe which bits should be considered first (continuing by
    # wrapping and going clockwise)
    interpolateDir = [
        [0b100011, 0b100111, 0b100111, 0b101111, 0b101111, 0b000000, 0b000000, 0b100011],
        [0b100011, 0b100111, 0b100111, 0b101111, 0b101111, 0b111111, 0b000000, 0b100011],
        [0b000000, 0b000100, 0b000100, 0b001100, 0b001100, 0b011100, 0b011100, 0b000000],
        [0b000000, 0b000100, 0b000100, 0b001100, 0b001100, 0b011100, 0b011100, 0b111100],
        [0b111001, 0b000000, 0b000000, 0b001001, 0b001001, 0b011001, 0b011001, 0b111001],
        [0b111001, 0b111101, 0b000000, 0b001001, 0b001001, 0b011001, 0b011001, 0b111001],
        [0b110010, 0b110110, 0b110110, 0b000000, 0b000000, 0b010010, 0b010010, 0b110010],
        [0b110010, 0b110110, 0b110110, 0b111110, 0b000000, 0b010010, 0b010010, 0b110010]
    ]

    # _surrounding_ will take the given indices and check the blob map
    # if any of the surrounding points adhere to the conditions stipulated by
    # the callable 'validator' argument. Two arrays  are returned, the first
    # array contains blobs that are edges (left, top, right, bottom), the second
    # array contains blobs that are corners (topleft, etc.).  Each value in the
    # returned lists contains the values [blobIndex, yIndex, xIndex]
    def _surrounding_(self, iY, iX, validator):
        # For readability the pixels are checked clockwise, starting at the
        # pixel to the left. The exception is the pixel common to all subcases
        # of the x-coordinate case
        blobsEdges = []
        blobsCorners = []

        if iX == 0:
            if validator(self.blobMap[iY, iX + 1]):
                blobsEdges.append([self.blobMap[iY, iX + 1], iY, iX + 1])

            if iY == 0:
                # Bottom-left corner
                if validator(self.blobMap[iY + 1, iX]):
                    blobsEdges.append([self.blobMap[iY + 1, iX], iY + 1, iX])

                if validator(self.blobMap[iY + 1, iX + 1]):
                    blobsCorners.append([self.blobMap[iY + 1, iX + 1], iY + 1, iX + 1])
            elif iY == len(self.axisY) - 1:
                # Top-left corner
                if validator(self.blobMap[iY - 1, iX + 1]):
                    blobsCorners.append([self.blobMap[iY - 1, iX + 1], iY - 1, iX + 1])

                if validator(self.blobMap[iY - 1, iX]):
                    blobsEdges.append([self.blobMap[iY - 1, iX], iY - 1, iX])
            else:
                # On the left border
                if validator(self.blobMap[iY + 1, iX]):
                    blobsEdges.append([self.blobMap[iY + 1, iX], iY + 1, iX])

                if validator(self.blobMap[iY + 1, iX + 1]):
                    blobsCorners.append([self.blobMap[iY + 1, iX + 1], iY + 1, iX])

                if validator(self.blobMap[iY - 1, iX + 1]):
                    blobsCorners.append([self.blobMap[iY - 1, iX + 1], iY - 1, iX + 1])

                if validator(self.blobMap[iY - 1, iX]):
                    blobsEdges.append([self.blobMap[iY - 1, iX], iY - 1, iX])
        elif iX == len(self.axisX) - 1:
            if validator(self.blobMap[iY, iX - 1]):
                blobsEdges.append([self.blobMap[iY, iX - 1], iY, iX - 1])

            if iY == 0:
                # bottom-right corner
                if validator(self.blobMap[iY + 1, iX - 1]):
                    blobsCorners.append([self.blobMap[iY + 1, iX - 1], iY + 1, iX - 1])

                if validator(self.blobMap[iY + 1, iX]):
                    blobsEdges.append([self.blobMap[iY + 1, iX], iY + 1, iX])
            elif iY == len(self.axisY) - 1:
                # top-right corner
                if validator(self.blobMap[iY - 1, iX]):
                    blobsEdges.append([self.blobMap[iY - 1, iX], iY - 1, iX])

                if validator(self.blobMap[iY - 1, iX - 1]):
                    blobsCorners.append([self.blobMap[iY - 1, iX - 1], iY - 1, iX - 1])
            else:
                # On the right border
                if validator(self.blobMap[iY + 1, iX - 1]):
                    blobsCorners.append([self.blobMap[iY + 1, iX - 1], iY + 1, iX - 1])

                if validator(self.blobMap[iY + 1, iX]):
                    blobsEdges.append([self.blobMap[iY + 1, iX], iY + 1, iX])

                if validator(self.blobMap[iY - 1, iX]):
                    blobsEdges.append([self.blobMap[iY - 1, iX], iY - 1, iX])

                if validator(self.blobMap[iY - 1, iX - 1]):
                    blobsCorners.append([self.blobMap[iY - 1, iX - 1], iY - 1, iX - 1])
        else:
            # X-coordinate is in the center
            if validator(self.blobMap[iY, iX - 1]):
                blobsEdges.append([self.blobMap[iY, iX - 1], iY, iX - 1])

            if validator(self.blobMap[iY, iX + 1]):
                blobsEdges.append([self.blobMap[iY, iX + 1], iY, iX + 1])

            if iY == 0:
                # On the bottom border
                if validator(self.blobMap[iY + 1, iX - 1]):
                    blobsCorners.append([self.blobMap[iY + 1, iX - 1], iY + 1, iX - 1])

                if validator(self.blobMap[iY + 1, iX]):
                    blobsEdges.append([self.blobMap[iY + 1, iX], iY + 1, iX])

                if validator(self.blobMap[iY + 1, iX + 1]):
                    blobsCorners.append([self.blobMap[iY + 1, iX + 1], iY + 1, iX + 1])
            elif iY == len(self.axisY) - 1:
                #On the top border
                if validator(self.blobMap[iY - 1, iX + 1]):
                    blobsCorners.append([self.blobMap[iY - 1, iX + 1], iY - 1, iX + 1])

                if validator(self.blobMap[iY - 1, iX]):
                    blobsEdges.append([self.blobMap[iY - 1, iX], iY - 1, iX])

                if validator(self.blobMap[iY - 1, iX - 1]):
                    blobsCorners.append([self.blobMap[iY - 1, iX - 1], iY - 1, iX - 1])
            else:
                # Somewhere in the center, all points can be added
                if validator(self.blobMap[iY + 1, iX - 1]):
                    blobsCorners.append([self.blobMap[iY + 1, iX - 1], iY + 1, iX - 1])

                if validator(self.blobMap[iY + 1, iX]):
                    blobsEdges.append([self.blobMap[iY + 1, iX], iY + 1, iX])

                if validator(self.blobMap[iY + 1, iX + 1]):
                    blobsCorners.append([self.blobMap[iY + 1, iX + 1], iY + 1, iX + 1])

                if validator(self.blobMap[iY - 1, iX + 1]):
                    blobsCorners.append([self.blobMap[iY - 1, iX + 1], iY - 1, iX + 1])

                if validator(self.blobMap[iY - 1, iX]):
                    blobsEdges.append([self.blobMap[iY - 1, iX], iY - 1, iX])

                if validator(self.blobMap[iY - 1, iX - 1]):
                    blobsCorners.append([self.blobMap[iY - 1, iX - 1], iY - 1, iX - 1])

        return blobsEdges, blobsCorners

    # _surroundingBlobs_ will check the blob map for existing blobs surrounding
    # the provided index
    def _surroundingBlobs_(self, iY, iX):
        return self._surrounding_(iY, iX, lambda val: val > 0)

    # _surroundingZeros_ will check the blob map for existing non-blobs
    # surrounding the provided index
    def _surroundingZeros_(self, iY, iX):
        return self._surrounding_(iY, iX, lambda val: val < 0)

    # _renumberBlobs_ will loop through all set blob indices. If the blob index
    # matches the original index passed as argument to this function then it
    # will be replaced by the new index
    def _renumberBlobs_(self, original, new):
        for iY in range(0, len(self.axisY)):
            for iX in range(0, len(self.axisX)):
                if self.blobMap[iY, iX] == original:
                    self.blobMap[iY, iX] = new

    # _nextEdge_ will take the given indices loop through all adjacent data
    # points (starting left, then continuing clockwise) to find a point that has
    # the same blob index. This function will return (yNew, xNew, direction) as
    # a tuple. The given direction corresponds to the constants defined at the
    # top of this file.
    def _nextEdge_(self, iY, iX, blobIndex, prevDirection):
        for i in range(prevDirection - 3, prevDirection + 5):
            iMod = i % 8
            checkY = iY + self.tileOffsets[iMod][0]
            checkX = iX + self.tileOffsets[iMod][1]

            if checkY < 0 or checkY >= len(self.axisY) or \
                    checkX < 0 or checkX >= len(self.axisX):
                continue

            if self.blobMap[checkY, checkX] == blobIndex:
                return checkY, checkX, iMod

        # Failed to find a surrounding data point, return this point
        return iY, iX, self.NONE

    # _interpolateGeneral_ contains the general point-interpolation code used
    # in _interpolateLeft_, _interpolateTop_, etc. without checking for
    # boundaries
    def _interpolateGeneral_(self, iY1, iX1, iY2, iX2, interpolant, selector):
        if self.data[iY2, iX2] < self.data[iY1, iX1] or self.blobMap[iY1, iX1] < 0:
            return TrackCommon.Lerp(self.data[iY1, iX1], interpolant[selector(iY1, iX1)],
                                    self.data[iY2, iX2], interpolant[selector(iY2, iX2)],
                                    self.vMin), self.blobMap[iY2, iX2]

        return TrackCommon.Lerp(self.data[iY1, iX1], interpolant[selector(iY1, iX1)],
                                self.data[iY2, iX2], interpolant[selector(iY2, iX2)],
                                self.vMax), self.blobMap[iY2, iX2]

    # _interpolateMultipleGeneral_ contains the general multi-datapoint
    # interpolation code used in _interpolateMultipleLeft_,
    # _interpolateMultipleTop_, etc. without checking for boundaries
    def _interpolateMultipleGeneral_(self, iY1, iX1, iY2, iX2, interpolant, selector):
        origin = interpolant[selector(iY1, iX1)]
        final = interpolant[selector(iY2, iX2)]
        delta = abs(origin - final)
        blobIndex = 0

        if self.blobMap[iY1, iX1] > 0:
            # The considered datapoint is part of a blob. This implies that the
            # current datapoints lies between all vMins and vMaxs for each of
            # the datasets
            for iData in range(0, self.numDatas):
                if self.datas[iData, iY2, iX2] < self.vMins[iData]:
                    new = TrackCommon.Lerp(self.datas[iData, iY1, iX1], interpolant[selector(iY1, iX1)],
                                           self.datas[iData, iY2, iX2], interpolant[selector(iY2, iX2)],
                                           self.vMins[iData])
                    newDelta = abs(origin - new)

                    if newDelta < delta:
                        final = new
                        delta = newDelta
                        blobIndex = self.blobMap[iY2, iX2]

                elif self.datas[iData, iY2, iX2] > self.vMaxs[iData]:
                    new = TrackCommon.Lerp(self.datas[iData, iY1, iX1], interpolant[selector(iY1, iX1)],
                                           self.datas[iData, iY2, iX2], interpolant[selector(iY2, iX2)],
                                           self.vMaxs[iData])
                    newDelta = abs(origin - new)

                    if newDelta < delta:
                        final = new
                        delta = newDelta
                        blobIndex = self.blobMap[iY2, iX2]

                # If the two cases above are not true then the current datamap
                # constituting the combined datamap is not providing an edge for
                # the combined datamap
        else:
            # The considered datpoint is part of a non-blob. This requires
            # checking where the datapoint is with respect to the vMin and the
            # vMax value
            for iData in range(0, self.numDatas):
                # TODO: Remove this check after I'm sure it is not necessary
                if self.datas[iData, iY2, iX2] < self.vMins[iData] and \
                        self.datas[iData, iY2, iX2] > self.vMaxs[iData]:
                    raise ValueError('Your logic is bad and you should feel bad. ' +
                        'Now put your brain to some proper use and figure out why ' +
                        'the non-blob is next to a non-blob')

                # Note that a non-blob datapoint might still reside in one of
                # the datamaps, but not in all datamaps
                if self.datas[iData, iY1, iX1] < self.vMins[iData]:
                    new = TrackCommon.Lerp(self.datas[iData, iY1, iX1], interpolant[selector(iY1, iX1)],
                                           self.datas[iData, iY2, iX2], interpolant[selector(iY2, iX2)],
                                           self.vMins[iData])
                    newDelta = abs(origin - new)

                    if newDelta < delta:
                        final = new
                        delta = newDelta
                        blobIndex = self.blobMap[iY2, iX2]

                elif self.datas[iData, iY1, iX1] > self.vMaxs[iData]:
                    new = TrackCommon.Lerp(self.datas[iData, iY1, iX1], interpolant[selector(iY1, iX1)],
                                           self.datas[iData, iY2, iX2], interpolant[selector(iY2, iX2)],
                                           self.vMaxs[iData])
                    newDelta = abs(origin - new)

                    if newDelta < delta:
                        final = new
                        delta = newDelta
                        blobIndex = self.blobMap[iY2, iX2]

        return final, blobIndex

    # _interpolateLeft_ will take the given indices and find the edge of the
    # contour on the left side
    def _interpolateLeft_(self, iY, iX):
        if iX == 0:
            return [self.axisX[iX], self.axisY[iY]], 0

        interpolated, blobIndex = self._interpolateGeneral_(iY, iX, iY, iX - 1, self.axisX, lambda y, x: x)
        return [ interpolated, self.axisY[iY] ], blobIndex

    # _interpolateTop_ will take the given indices and find the edge of the
    # contour on the top side
    def _interpolateTop_(self, iY, iX):
        if iY == len(self.axisY) - 1:
            return [self.axisX[iX], self.axisY[iY]], 0

        interpolated, blobIndex = self._interpolateGeneral_(iY, iX, iY + 1, iX, self.axisY, lambda y, x: y)
        return [ self.axisX[iX], interpolated ], blobIndex

    # _interpolateRight_ will take the given indices and find the edge of the
    # contour on the right side
    def _interpolateRight_(self, iY, iX):
        if iX == len(self.axisX) - 1:
            return [self.axisX[iX], self.axisY[iY]], 0

        interpolated, blobIndex = self._interpolateGeneral_(iY, iX, iY, iX + 1, self.axisX, lambda y, x: x)
        return [ interpolated, self.axisY[iY] ], blobIndex

    # _interpolateBottom_ will take the given indices and find the edge of the
    # contour on the bottom side
    def _interpolateBottom_(self, iY, iX):
        if iY == 0:
            return [self.axisX[iX], self.axisY[iY]], 0

        interpolated, blobIndex = self._interpolateGeneral_(iY, iX, iY - 1, iX, self.axisY, lambda y, x: y)
        return [ self.axisX[iX], interpolated ], blobIndex

    # _interpolateMultipleLeft_ will take the given indices and find the edge of
    # the contour on the left side, taking the multiple maps into account by
    # finding the edge closest to the center of the supplied indices of all maps
    # constituting to the edge
    def _interpolateMultipleLeft_(self, iY, iX):
        if iX == 0:
            return [self.axisX[iX], self.axisY[iY]], 0

        interpolated, blobIndex = self._interpolateMultipleGeneral_(iY, iX, iY, iX - 1, self.axisX, lambda y, x: x)
        return [ interpolated, self.axisY[iY] ], blobIndex

    # _interpolateMultipleTop_ will take the given indices and find the edge of
    # the contour on the top side, taking multiple maps into account by finding
    # the edge closest to the center of the supplied indices of all maps
    # constituting the edge
    def _interpolateMultipleTop_(self, iY, iX):
        if iY == len(self.axisY) - 1:
            return [self.axisX[iX], self.axisY[iY]], 0

        interpolated, blobIndex = self._interpolateMultipleGeneral_(iY, iX, iY + 1, iX, self.axisY, lambda y, x: y)
        return [ self.axisX[iX], interpolated ], blobIndex

    # _interpolateMultipleRight_ will take the given indices and find the edge
    # of the contour on the right side, taking multiple maps into account by
    # finding the edge closest to the center of the supplied indices of all maps
    # constituting the edge
    def _interpolateMultipleRight_(self, iY, iX):
        if iX == len(self.axisX) - 1:
            return [self.axisX[iX], self.axisY[iY]], 0

        interpolated, blobIndex = self._interpolateMultipleGeneral_(iY, iX, iY, iX + 1, self.axisX, lambda y, x: x)
        return [ interpolated, self.axisY[iY] ], blobIndex

    # _interpolateMultipleBottom_ will take the given indices and find the edge
    # of the contour on the bottom side, taking multiple maps into account by
    # finding the edge closest to the center of the supplied indices of all maps
    # constituting the edge
    def _interpolateMultipleBottom_(self, iY, iX):
        if iY == 0:
            return [self.axisX[iX], self.axisY[iY]], 0

        interpolated, blobIndex = self._interpolateMultipleGeneral_(iY, iX, iY - 1, iX, self.axisY, lambda y, x: y)
        return [ self.axisX[iX], interpolated ], blobIndex

    def _generateContours_(self):
        self.blobMap = np.zeros(self.data.shape, dtype=np.int16)

        # blobSets will keep track of a single edge pixel belonging to a single
        # blob. Each entry in the blobSets list contains
        # [blobIndex, yIndex, xIndex, [yZeroIndex, xZeroIndex]]
        blobCounter = 0
        zeroCounter = 0
        blobSets = []

        for iY in range(0, len(self.axisY)):
            for iX in range(0, len(self.axisX)):
                blobsEdges = None
                blobsCorners = None

                # Check if any existing blobs surround the current pixels.
                # From those determine the blob number, create a new blob
                # number and optionally renumber an existing blob
                isBlob = self.data[iY,iX] >= self.vMin and \
                    self.data[iY,iX] <= self.vMax

                if isBlob:
                    blobsEdges, blobsCorners = self._surroundingBlobs_(iY, iX)
                else:
                    blobsEdges, blobsCorners = self._surroundingZeros_(iY, iX)

                chosenBlob = 0
                reblobbed = []

                for i in range(0, len(blobsEdges)):
                    if chosenBlob == 0:
                        # Found a blob number
                        chosenBlob = blobsEdges[i][0]
                    elif blobsEdges[i][0] != chosenBlob and not (blobsEdges[i][0] in reblobbed):
                        # Found another adjacent blob whose number doesn't
                        # match the assigned one: renumber the blob
                        self._renumberBlobs_(blobsEdges[i][0], chosenBlob)
                        reblobbed.append(blobsEdges[i][0])

                        # Adjust the number in the blobSets array
                        for iBlob in range(0, len(blobSets)):
                            if blobSets[iBlob][0] == blobsEdges[i][0]:
                                #blobSets[iBlob][0] = chosenBlob
                                del blobSets[iBlob]
                                break

                for i in range(0, len(blobsCorners)):
                    # Check if the corner blob can be considered part of the
                    # current blob by checking the value of the point
                    # inbetween four corners
                    midValue = (
                        self.data[iY, iX] +
                        self.data[iY, blobsCorners[i][2]] +
                        self.data[blobsCorners[i][1], iX] +
                        self.data[blobsCorners[i][1], blobsCorners[i][2]]
                    ) / 4

                    if (midValue >= self.vMin and midValue <= self.vMax) or (not isBlob):
                        # Consider the blob to be part of the currently
                        # considered data value
                        if chosenBlob == 0:
                            chosenBlob = blobsCorners[i][0]
                        elif blobsCorners[i][0] != chosenBlob and not (blobsCorners[i][0] in reblobbed):
                            self._renumberBlobs_(blobsCorners[i][0], chosenBlob)
                            reblobbed.append(blobsCorners[i][0])

                            for iBlob in range(0, len(blobSets)):
                                if blobSets[iBlob][0] == blobsCorners[i][0]:
                                    #blobSets[iBlob][0] = chosenBlob
                                    del blobSets[iBlob]
                                    break

                if chosenBlob == 0:
                    # The current data point does not belong to another blob
                    if isBlob:
                        blobCounter += 1
                        chosenBlob = blobCounter
                        blobSets.append([chosenBlob, iY, iX])
                    else:
                        zeroCounter -= 1
                        chosenBlob = zeroCounter
                        blobSets.append([chosenBlob, iY, iX])

                self.blobMap[iY, iX] = chosenBlob

        # Use the generated blobSets variable to find the borders of each blob.
        directionMap = np.empty(self.data.shape, dtype=np.uint8)
        directionMap.fill(self.NONE)

        includers = []
        excluders = []

        for iBlob in range(0, len(blobSets)):
            # Perform the first two iterations to get initial values for the
            # continuous algorithm
            curContour = []
            curBlobIndex = blobSets[iBlob][0]
            targetBlob = 0

            curY, curX, curDirection = self._nextEdge_(blobSets[iBlob][1],
                blobSets[iBlob][2], curBlobIndex, self.START)
            newY, newX, newDirection = self._nextEdge_(curY, curX, curBlobIndex, curDirection)

            # If the blob contains values outside of the indicated bound then
            # it should be disregarded if any of the points is on the edge
            isNonBlobAndEdge = False

            if curBlobIndex < 0 and (
                    curX == 0 or curX == len(self.axisX) - 1 or
                    curY == 0 or curY == len(self.axisY) - 1 or
                    newX == 0 or newX == len(self.axisX) - 1 or
                    newY == 0 or newY == len(self.axisY) - 1):
                continue

            # Special case: A single value
            if curX == newX and curY == newY:
                point, blobIndex = self._interpolateLeft_(curY, curX)
                curContour.append(point)

                if blobIndex != 0:
                    targetBlob = blobIndex

                point, blobIndex = self._interpolateTop_(curY, curX)
                curContour.append(point)

                if blobIndex != 0:
                    targetBlob = blobIndex

                point, blobIndex = self._interpolateRight_(curY, curX)
                curContour.append(point)

                if blobIndex != 0:
                    targetBlob = blobIndex

                point, blobIndex = self._interpolateBottom_(curY, curX)
                curContour.append(point)

                if blobIndex != 0:
                    targetBlob = blobIndex
            else:
                # Loop through all other points
                while directionMap[newY, newX] != newDirection:
                    if directionMap[newY, newX] == self.NONE:
                        # First time entering this data point, set the direction in
                        # which it was entered
                        directionMap[newY, newX] = newDirection

                    # Use the direction map to figure out which interpolations to
                    # perform
                    interpolate = self.interpolateDir[curDirection][newDirection]

                    dirs = (interpolate >> 2) & 0b1111
                    shift = interpolate & 0b11
                    mask = 0b0001 << shift

                    for i in range(0, 4):
                        if dirs & mask:
                            point = None
                            blobIndex = 0

                            if mask == self.BITMAP_LEFT:
                                point, blobIndex = self._interpolateLeft_(curY, curX)
                            elif mask == self.BITMAP_TOP:
                                point, blobIndex = self._interpolateTop_(curY, curX)
                            elif mask == self.BITMAP_RIGHT:
                                point, blobIndex = self._interpolateRight_(curY, curX)
                            elif mask == self.BITMAP_BOTTOM:
                                point, blobIndex = self._interpolateBottom_(curY, curX)
                            else:
                                # TODO: Remove this if it is certain that this
                                # piece of logic is okay
                                raise ValueError("You suck at considering cases, " +
                                    "Max. You really suck at it... Hard! You should " +
                                    "be ashamed of your existence.")

                            curContour.append(point)

                            if blobIndex != 0:
                                targetBlob = blobIndex

                        mask <<= 1

                        if mask == 0b10000:
                            mask = 0b0001

                    # Shift searching variables for next iteration
                    curX = newX
                    curY = newY
                    curDirection = newDirection

                    # Check if this is a non-blob on the edge
                    newY, newX, newDirection = self._nextEdge_(curY, curX, curBlobIndex, curDirection)

                    if curBlobIndex < 0 and (
                            curX == 0 or curX == len(self.axisX) - 1 or
                            curY == 0 or curY == len(self.axisY) - 1):
                        isNonBlobAndEdge = True
                        break

            # Add the newly found contour to the array of contour after closing
            # it using the initially found point
            if not isNonBlobAndEdge:
                curContour.append(curContour[0]) # close the contour

                if curBlobIndex < 0:
                    excluders.append([targetBlob, np.asarray(curContour)])
                else:
                    includers.append([curBlobIndex, np.asarray(curContour)])

        # Construct includers with the corresponding excluders
        for iIncluder in range(0, len(includers)):
            newContour = ContourSingleSection(self,
                                              includers[iIncluder][1],
                                              includers[iIncluder][0])

            for iExcluder in range(0, len(excluders)):
                if excluders[iExcluder][0] == includers[iIncluder][0]:
                    newContour.addExcluder(excluders[iExcluder][1],
                                           excluders[iExcluder][0])

            self.contours.append(newContour)

    def _combineContours_(self):
        # AND all the data maps together into a single blobmap
        self.blobMap = np.zeros([self.datas.shape[1], self.datas.shape[2]], dtype=np.int16)

        blobCounter = 0
        zeroCounter = 0
        blobSets = []

        for iY in range(0, len(self.axisY)):
            for iX in range(0, len(self.axisX)):
                # Find surrounding blobs (or non-blobs)
                blobsEdges = None
                blobsCorners = None

                isBlob = True

                for iData in range(0, self.numDatas):
                    if self.datas[iData, iY, iX] < self.vMins[iData] or \
                            self.datas[iData, iY, iX] > self.vMaxs[iData]:
                        isBlob = False

                if isBlob:
                    blobsEdges, blobsCorners = self._surroundingBlobs_(iY, iX)
                else:
                    blobsEdges, blobsCorners = self._surroundingZeros_(iY, iX)

                chosenBlob = 0
                reblobbed = []

                # Loop through the surrounding blobs to look for a pre-existing
                # blob number (or zero-number)
                for i in range(0, len(blobsEdges)):
                    if chosenBlob == 0:
                        chosenBlob = blobsEdges[i][0]
                    elif blobsEdges[i][0] != chosenBlob and not (blobsEdges[i][0] in reblobbed):
                        # Renumber the previous blob number
                        self._renumberBlobs_(blobsEdges[i][0], chosenBlob)
                        reblobbed.append(blobsEdges[i][0])

                        for iBlob in range(0, len(blobSets)):
                            if blobSets[iBlob][0] == blobsEdges[i][0]:
                                del blobSets[iBlob]
                                break

                for i in range(0, len(blobsCorners)):
                    # Check all data maps to see if the center point belongs
                    # to the current blob
                    iXAlt = blobsCorners[i][2]
                    iYAlt = blobsCorners[i][1]

                    isPartOfBlob = True

                    for iData in range(0, self.numDatas):
                        midValue = (
                            self.datas[iData, iY, iX] +
                            self.datas[iData, iY, iXAlt] +
                            self.datas[iData, iYAlt, iX] +
                            self.datas[iData, iYAlt, iXAlt]
                        ) / 4

                        if midValue < self.vMins[iData] or midValue > self.vMaxs[iData]:
                            isPartOfBlob = False
                            break

                    if isPartOfBlob or (not isBlob):
                        # Consider the considered point to be part of the
                        # current blob
                        if chosenBlob == 0:
                            chosenBlob = blobsCorners[i][0]
                        elif blobsCorners[i][0] != chosenBlob and not (blobsCorners[i][0] in reblobbed):
                            self._renumberBlobs_(blobsCorners[i][0], chosenBlob)
                            reblobbed.append(blobsCorners[i][0])

                            for iBlob in range(0, len(blobSets)):
                                if blobSets[iBlob][0] == blobsCorners[i][0]:
                                    del blobSets[iBlob]
                                    break

                if chosenBlob == 0:
                    # Did not find a blob to which this point belongs
                    if isBlob:
                        blobCounter += 1
                        chosenBlob = blobCounter
                        blobSets.append([chosenBlob, iY, iX])
                    else:
                        zeroCounter -= 1
                        chosenBlob = zeroCounter
                        blobSets.append([chosenBlob, iY, iX])

                self.blobMap[iY, iX] = chosenBlob

        # Use the generate blobSets variable to find the borders of the
        # compound blob. All underlying datasets should be interpolated if
        # they're part of the compound blob's edges
        directionMap = np.empty([self.datas.shape[1], self.datas.shape[2]], dtype=np.uint8)
        directionMap.fill(self.NONE)

        includers = []
        excluders = []

        for iBlob in range(0, len(blobSets)):
            # Perform the first two iterations to get initial values for the
            # continous algorithm
            curContour = []
            curBlobIndex = blobSets[iBlob][0]
            targetBlob = 0

            curY, curX, curDirection = self._nextEdge_(blobSets[iBlob][1],
                blobSets[iBlob][2], curBlobIndex, self.START)

            newY, newX, newDirection = self._nextEdge_(curY, curX, curBlobIndex, curDirection)

            # If a non-blob contains values that are on the border of the data
            # map then it should be disregarded
            isNonBlobAndEdge = False

            if curBlobIndex < 0 and (
                    curX == 0 or curX == len(self.axisX) - 1 or
                    curY == 0 or curY == len(self.axisY) - 1 or
                    newX == 0 or newX == len(self.axisX) - 1 or
                    newY == 0 or newY == len(self.axisY) - 1):
                continue

            # Special case: A single value
            if curX == newX and curY == newY:
                point, blobIndex = self._interpolateMultipleLeft_(curY, curX)
                curContour.append(point)

                if blobIndex != 0:
                    targetBlob = blobIndex

                point, blobIndex = self._interpolateMultipleTop_(curY, curX)
                curContour.append(point)

                if blobIndex != 0:
                    targetBlob = blobIndex

                point, blobIndex = self._interpolateMultipleRight_(curY, curX)
                curContour.append(point)

                if blobIndex != 0:
                    targetBlob = blobIndex

                point, blobIndex = self._interpolateMultipleBottom_(curY, curX)
                curContour.append(point)

                if blobIndex != 0:
                    targetBlob = blobIndex
            else:
                # Walk around the contour
                while directionMap[newY, newX] != newDirection:
                    if directionMap[newY, newX] == self.NONE:
                        # First time entering this point, set the direction in
                        # which it was entered
                        directionMap[newY, newX] = newDirection

                    # Perform interpolation on the edge
                    interpolate = self.interpolateDir[curDirection][newDirection]

                    dirs = (interpolate >> 2) & 0b1111
                    shift = interpolate & 0b11
                    mask = 0b0001 << shift

                    for i in range(0, 4):
                        if dirs & mask:
                            point = None
                            blobIndex = 0

                            if mask == self.BITMAP_LEFT:
                                point, blobIndex = self._interpolateMultipleLeft_(curY, curX)
                            elif mask == self.BITMAP_TOP:
                                point, blobIndex = self._interpolateMultipleTop_(curY, curX)
                            elif mask == self.BITMAP_RIGHT:
                                point, blobIndex = self._interpolateMultipleRight_(curY, curX)
                            elif mask == self.BITMAP_BOTTOM:
                                point, blobIndex = self._interpolateMultipleBottom_(curY, curX)
                            else:
                                # TODO: Remove this when it is sure that you don't need it
                                raise ValueError("Another check to ensure you're not " +
                                    "being an idiot. Just so you know: Your past self " +
                                    "assumed you were an idiot. That's probably why " +
                                    "you see this message. If you think yourself now, " +
                                    "then you will've thanked your past self in the " +
                                    "future.")

                            curContour.append(point)

                            if blobIndex != 0:
                                targetBlob = blobIndex

                        mask <<= 1

                        if mask == 0b10000:
                            mask = 0b0001

                    # Prepare variables for the next edge point iteration
                    curX = newX
                    curY = newY
                    curDirection = newDirection

                    # Check if this is a non-blob on the edge
                    newY, newX, newDirection = self._nextEdge_(curY, curX, curBlobIndex, curDirection)

                    if curBlobIndex < 0 and (
                            curX == 0 or curX == len(self.axisX) - 1 or
                            curY == 0 or curY == len(self.axisY) - 1):
                        isNonBlobAndEdge = True
                        break

            if not isNonBlobAndEdge:
                curContour.append(curContour[0]) # close the contour

                if curBlobIndex < 0:
                    excluders.append([targetBlob, np.asarray(curContour)])
                else:
                    includers.append([curBlobIndex, np.asarray(curContour)])

        # Construct includers with the corresponding excluders
        for iIncluder in range(0, len(includers)):
            newContour = ContourSingleSection(self,
                                              includers[iIncluder][1],
                                              includers[iIncluder][0])

            for iExcluder in range(0, len(excluders)):
                if excluders[iExcluder][0] == includers[iIncluder][0]:
                    newContour.addExcluder(excluders[iExcluder][1],
                                           excluders[iExcluder][0])

            self.contours.append(newContour)

    def _reset_(self):
        self.axisX = []
        self.axisXLookup = None
        self.axisY = []
        self.axisYLookup = None
        self.data = []
        self.vMin = 0.0
        self.vMax = 0.0
        self.contours = []
        self.blobMap = None
        self.edgeMap = None

    def _setAxes_(self, axisX, axisY):
        if TrackCommon.IsAscending(axisX):
            self.axisXLookup = TrackCommon.Find1DBisectionAscending
        elif TrackCommon.IsDescending(axisX):
            self.axisXLookup = TrackCommon.Find1DBisectionDescending
        else:
            raise ValueError("axisX is neither ascending nor descending")

        self.axisX = np.asarray(axisX)

        if TrackCommon.IsAscending(axisY):
            self.axisYLookup = TrackCommon.Find1DBisectionAscending
        elif TrackCommon.IsDescending(axisY):
            self.axisYLookup = TrackCommon.Find1DBisectionDescending
        else:
            raise ValueError("axisY is neither ascending nor descending")

        self.axisY = np.asarray(axisY)

    def __init__(self):
        self._reset_()

    def setData(self, axisX, axisY, data, vMin, vMax):
        self._reset_()

        # Preshape and check the input data
        if not TrackCommon.IsArray(axisX):
            raise ValueError("Expected 'axisX' to be an array")

        if not TrackCommon.IsArray(axisY):
            raise ValueError("Expected 'axisY' to be an array")

        if not TrackCommon.IsArray(data):
            raise ValueError("Expected 'data' to be an array")

        self._setAxes_(axisX, axisY)
        data = np.asarray(data)

        if len(axisX.shape) != 1:
            raise ValueError("Expected 'axisX' to be a 1D array")

        if len(axisY.shape) != 1:
            raise ValueError("Expected 'axisY' to be a 1D array")

        if len(data.shape) != 2:
            raise ValueError("Expected data to be a 2D array")

        if data.shape[0] != axisY.shape[0]:
            raise ValueError("Expected first data axis to correspond to y-axis length")

        if data.shape[1] != axisX.shape[0]:
            raise ValueError("Expected second data axis to correspond to x-axis length")

        if vMin > vMax:
            raise ValueError("Expected 'vMin' to be smaller than 'vMax'")

        # Store values
        self.axisX = axisX
        self.axisY = axisY
        self.data = data
        self.vMin = vMin
        self.vMax = vMax

        # Perform contour-finding algorithm
        self._generateContours_()

    def combineData(self, axisX, axisY, datas, vMins, vMaxs):
        self._reset_()

        # Preshape and check the input data
        if not TrackCommon.IsArray(axisX):
            raise ValueError("Expected 'axisX' to be an array")

        if not TrackCommon.IsArray(axisY):
            raise ValueError("Expected 'axisY' to be an array")

        if not TrackCommon.IsArray(vMins):
            raise ValueError("Expected 'vMins' to be an array")

        if not TrackCommon.IsArray(vMaxs):
            raise ValueError("Expected 'vMaxs' to be an array")

        self._setAxes_(axisX, axisY)
        self.vMins = np.asarray(vMins)
        self.vMaxs = np.asarray(vMaxs)

        if len(self.axisX.shape) != 1:
            raise ValueError("Expected 'axisX' to be a 1D array")

        if len(self.axisY.shape) != 1:
            raise ValueError("Expected 'axisY' to be a 1D array")

        if len(self.vMins.shape) != 1:
            raise ValueError("Expected 'vMins' to be a 1D array")

        if len(self.vMaxs.shape) != 1:
            raise ValueError("Expected 'vMaxs' to be a 1D array")

        self.datas = np.asarray(datas)
        self.numDatas = len(self.datas)

        if self.datas.shape[0] != self.vMins.shape[0]:
            raise ValueError("Expected 'datas' first dimension length to equal 'vMins' length")

        if self.datas.shape[0] != self.vMaxs.shape[0]:
            raise ValueError("Expected 'datas' first dimension length to equal 'vMaxs' length")

        if self.datas.shape[1] != self.axisY.shape[0]:
            raise ValueError("Expected 'datas' second dimension length to equal 'axisY' length")

        if self.datas.shape[2] != self.axisX.shape[0]:
            raise ValueError("Expected 'datas' third dimension length to equal 'axisX' length")

        self._combineContours_()

    def getNumContours(self):
        return len(self.contours)

    def getContour(self, index):
        return self.contours[index]

import matplotlib.pyplot as plt
import matplotlib as mpl

def __TestContourPlot__(contour):
    # Retrieve minimum and maximum values
    fig = plt.figure()
    xMin = np.min(contour.axisX)
    xMax = np.max(contour.axisX)
    yMin = np.min(contour.axisY)
    yMax = np.max(contour.axisY)

    vMin = None
    vMax = None

    if len(contour.data) == 0:
        vMin = np.min(contour.datas)
        vMax = np.max(contour.datas)
    else:
        vMin = np.min(contour.data)
        vMax = np.max(contour.data)

    norm = mpl.colors.Normalize(vMin, vMax)
    deltaX = (contour.axisX[1] - contour.axisX[0]) / 2.0
    deltaY = (contour.axisY[1] - contour.axisY[0]) / 2.0

    ax = fig.add_subplot(121)
    cmap = plt.get_cmap('jet')

    # Plot datamap
    if len(contour.data) == 0:
        for i in range(0, contour.numDatas):
            ax.imshow(contour.datas[i], extent=[xMin - deltaX, xMax + deltaX,
                yMin - deltaY, yMax + deltaY], norm=norm, cmap=cmap, origin='lower',
                alpha=(1.0 / contour.numDatas))
    else:
        ax.imshow(contour.data, extent=[xMin - deltaX, xMax + deltaX,
            yMin - deltaY, yMax + deltaY], norm=norm, cmap=cmap, origin='lower')
    ax.grid(True)

    # Plot contours
    for i in range(0, len(contour.contours)):
        curContour = contour.contours[i]
        includer = curContour.getIncluder()
        ax.plot(includer[:, 0], includer[:, 1], 'g')

        for j in range(0, curContour.getNumExcluders()):
            excluder = curContour.getExcluder(j)
            ax.plot(excluder[:, 0], excluder[:, 1], 'r')

    ax = fig.add_subplot(122)
    ax.imshow(contour.blobMap, extent=[xMin - deltaX, xMax + deltaX,
        yMin - deltaY, yMax + deltaY], norm=norm, cmap=cmap, origin='lower')

    # Plot horizontal ranges
    for i in range(1, len(contour.axisY)):
        curY = (contour.axisY[i - 1] + contour.axisY[i]) / 2.0

        for j in range(0, len(contour.contours)):
            curContour = contour.contours[j]
            ranges = curContour.getHorizontalRanges(curY)

            for k in range(0, len(ranges)):
                ax.plot([ranges[k, 0], ranges[k, 1]], [curY, curY], 'w')

    # Plot vertical ranges
    for i in range(1, len(contour.axisX)):
        curX = (contour.axisX[i - 1] + contour.axisX[i]) / 2.0

        for j in range(0, len(contour.contours)):
            curContour = contour.contours[j]
            ranges = curContour.getVerticalRanges(curX)

            for k in range(0, len(ranges)):
                ax.plot([curX, curX], [ranges[k, 0], ranges[k, 1]], 'w')

    ax.grid(True)

def __TestContourSimple__():
    contour = Contour()
    x = np.linspace(0, 4, 5)
    y = np.linspace(5, 8, 4)
    d = [
        [0, 1, 2, 1, 1],
        [1, 2, 4, 4, 1],
        [1, 4, 5, 4, 3],
        [1, 3, 0, 1, 2]
    ]

    contour.setData(x, y, d, 3.5, 6)
    __TestContourPlot__(contour)

def __TestContourHoles__():
    contour = Contour()
    x = np.linspace(0, 7, 8)
    y = np.linspace(0, 7, 8)
    d = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 1, 1, 0],
        [0, 1, 0, 0, 0, 0, 1, 0],
        [0, 1, 0, 1, 1, 0, 1, 0],
        [0, 1, 0, 1, 1, 0, 1, 0],
        [0, 1, 0, 0, 0, 0, 1, 0],
        [0, 1, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]

    contour.setData(x, y, d, 0.5, 2)
    __TestContourPlot__(contour)

def __TestContourEdges__():
    contour = Contour()

    x = np.linspace(0, 10, 11)
    y = np.linspace(0, 8, 9)
    d = [
        [0, 0, 0, 1, 2, 0, 2, 1, 0, 0, 0],
        [0, 0, 0, 2, 3, 0, 3, 2, 0, 0, 0],
        [2, 3, 0, 0, 1, 3, 2, 0, 0, 1, 3],
        [1, 4, 1, 0, 0, 0, 0, 0, 3, 2, 4],
        [0, 0, 2, 0, 0, 0, 0, 0, 3, 0, 0],
        [1, 4, 1, 0, 0, 0, 0, 0, 3, 2, 4],
        [2, 3, 0, 0, 1, 3, 2, 0, 0, 1, 3],
        [0, 0, 0, 2, 3, 0, 3, 2, 0, 0, 0],
        [0, 0, 0, 1, 2, 0, 2, 1, 0, 0, 0],
    ]

    contour.setData(x, y, d, 0.5, 6)
    __TestContourPlot__(contour)

def __TestContourNested__():
    contour = Contour()

    x = np.linspace(0, 10, 11)
    y = np.linspace(0, 10, 11)
    d = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0],
        [0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0],
        [0, 1, 0, 1, 0, 3, 0, 1, 0, 1, 0],
        [0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0],
        [0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ]

    for iY in range(0, 11):
        for iX in range(0, 11):
            if d[iY][iX] == 0:
                d[iY][iX] += 1.0 * np.random.rand() - 0.55

    contour.setData(x, y, d, 0.5, 2)
    __TestContourPlot__(contour)

def __TestContourMultipleNested__():
    contour = Contour()

    x = np.linspace(0, 6, 7)
    y = np.linspace(0, 6, 7)
    d = [
        [0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 1, 1],
        [0, 1, 0, 1, 0, 1, 0],
        [0, 1, 1, 1, 1, 1, 1],
        [0, 1, 0, 1, 0, 1, 0],
        [0, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0]
    ]

    contour.setData(x, y, d, 0.5, 1.5)
    __TestContourPlot__(contour)

def __TestContourMinMax__():
    contour = Contour()

    x = np.linspace(0, 5, 6)
    y = np.linspace(0, 5, 6)
    d = [
        [0, 1, 4, 3, 2, 5],
        [0, 2, 5, 2, 5, 3],
        [0, 1, 3, 4, 4, 1],
        [0, 3, 4, 5, 5, 4],
        [0, 6, 7, 7, 6, 8],
        [9, 10, 11, 8, 8, 10]
    ]

    contour.setData(x, y, d, 3, 5.5)
    __TestContourPlot__(contour)

def __TestContourCombineSimple__():
    contour = Contour()

    x = np.linspace(0, 5, 6)
    y = np.linspace(0, 5, 6)
    d1 = [
        [1, 1, 0, 0, 0, 0],
        [1, 1, 1, 0, 0, 0],
        [1, 1, 1, 1, 0, 0],
        [1, 1, 1, 1, 0, 0],
        [1, 1, 1, 0, 0, 0],
        [1, 1, 0, 0, 0, 0]
    ]
    d2 = [
        [0, 0, 0, 1, 1, 1],
        [0, 0, 1, 1, 1, 1],
        [0, 1, 1, 1, 1, 1],
        [0, 1, 1, 1, 1, 1],
        [0, 0, 1, 1, 1, 1],
        [0, 0, 0, 1, 1, 1]
    ]

    contour.combineData(x, y, [d1, d2], [0.5, 0.5], [1.5, 1.5])
    __TestContourPlot__(contour)

__TestContourSimple__()
__TestContourHoles__()
__TestContourEdges__()
__TestContourNested__()
__TestContourMultipleNested__()
#__TestContourMinMax__()
#__TestContourCombineSimple__()
