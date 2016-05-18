# -*- coding: utf-8 -*-
"""
Created on Wed May 18 22:54:20 2016

@author: MaxHenger
"""

def findClosestUnordered(array, toFind):
    delta = abs(toFind - array[0])
    best = 0

    for i in range(1, len(array)):
        curDelta = abs(toFind - array[i])

        if curDelta < delta:
            best = i
            delta = curDelta

    return best

def findClosestAscending(array, toFind):
    iLeft = 0
    iRight = len(array) - 1

    while True:
        iCenter = int((iLeft + iRight) / 2)
        current = array[iCenter]

        if current <= toFind:
            if iCenter - iLeft <= 1:
                if abs(toFind - array[iLeft]) < iCenter:
                    return iLeft

                return iCenter

            iRight = iCenter
        else:
            if iRight - iCenter <= 1:
                if abs(toFind - array[iRight]) < iCenter:
                    return iRight

                return iCenter

            iLeft = iCenter

def stripSpaces(string):
    result = ''

    for i in range(0, len(string)):
        if string[i] != ' ':
            result += string[i]

    return result
