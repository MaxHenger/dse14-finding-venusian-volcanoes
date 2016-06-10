# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 21:20:11 2016

@author: MaxHenger
"""

import time

class TimeEstimator:
    def __init__(self, totalIterations):
        self.start = 0
        self.iterationStart = 0
        self.stop = 0
        self.total = totalIterations
        self.current = 0

    def __padFront__(self, string, character, total):
        string = str(string)
        return character * (total - len(string)) + string

    def __formatTicks__(self, delta):
        delta = int(delta * 1e3)
        msecs = delta % 1000
        secs = (delta - msecs) / 1000
        minutes = secs / 60
        hours = minutes / 60
        secs = int(secs) % 60
        minutes = int(minutes) % 60
        hours = int(hours) % 24

        return self.__padFront__(hours, '0', 2) + ":" + \
            self.__padFront__(minutes, '0', 2) + ":" + \
            self.__padFront__(secs, '0', 2) + "." + \
            self.__padFront__(msecs, '0', 3)

    def startTiming(self):
        self.start = time.clock()

    def startIteration(self, iteration):
        self.iterationStart = time.clock()

    def finishedIteration(self, iteration):
        self.current = iteration
        self.stop = time.clock()

    def getTotalElapsed(self):
        return self.__formatTicks__(self.stop - self.start)

    def getIterationElapsed(self):
        return self.__formatTicks__(self.stop - self.iterationStart)

    def getEstimatedRemaining(self):
        spent = self.stop - self.start
        estimation = spent / (self.current + 1) * (self.total - 1 - self.current)
        return self.__formatTicks__(estimation)

    def getEstimatedTotal(self):
        spent = self.stop - self.start
        estimation = spent / (self.current + 1) * (self.total)
        return self.__formatTicks__(estimation)

    def getEstimatedEnd(self):
        spent = self.stop - self.start
        estimation = spent / (self.current + 1) * (self.total - 1 - self.current)
        end = time.mktime(time.localtime()) + estimation
        return self.__formatTicks__(end)
