# -*- coding: utf-8 -*-
"""
Created on Tue May  3 12:29:47 2016

This file contains the standard definition of an 'Atmosphere' as a class with
methods to obtain pressure, density, temperature, and zonal and meridional
wind speeds. Each method includes a way to retrieve uncertainties (in the form
of in minimum, mean and a maximum).

As general advise: always include the height, latitude and solar longitude if
they can be calculated. At the moment the solar longitude is not used at all
and the latitude is only used for the upper atmosphere. However, if your
written program is using these values correctly then an update to this file to
include a better model will immediately improve your written program.

In general all functions accept and return values in SI units. An exception is
angles, these are entered in degrees. Also attempt to call the functions with
arrays as much as possible to reduce unneccesary calculations.

The general use of this file is as following:
    
    include Atmosphere
    atm = Atmosphere.Atmosphere('preliminary')
    minPressure, meanPressure, maxPressure = atm.pressure(11000, 10, 10)
    
@author: MaxHenger
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as scp_ip
import Utility as util
import AtmosphereConstants as atm_const

# __determineInterpolationConstants1D__ is a helper functions to pre-calculate 
# 1D interpolation coefficients required for interpolation between the 
# atmospheric properties of Venus. This way the interpolation algorithm 
# deriving these coefficients does not have to be run each time a certain value 
# is required.
#
# Arguments
#   z: height (or: interpolation x-value)
#   v: parameter to interpolate (or: interpolation y-value)
#
# Outputs
#   knots: the knots of the spline
#   coeffs: the coefficients of the spline
#   order: the order (or: degree) of the spline
def __determineInterpolationConstants1D__(z, v, degree=None, ax=None):
    # Construct spline
    spline = []
    
    if degree is None:
        degree = 3
    
    spline = scp_ip.splrep(z, v, k=degree)
        
    # Draw spline if so required
    if not (ax is None):
        zRefined = np.linspace(z[0], z[-1], len(z) * 100)
        vRefined = scp_ip.splev(zRefined, spline)
        ax.plot(zRefined, vRefined)
     
    # Return coefficients
    return spline

# __determineInterpolationConstants2D__ is a hlper function to pre-calculate 2D
# interpolation coefficients required for interpolation between the atmospheric
# properties of Venus. This way the interpolation algorithm deriving these 
# coefficients does not have to be run each time a certain value is required
#
# Arguments:
#   lat: latitude (or: interpolation x-value)
#   z: height (or: interpolation y-value)
#   v: parameter to interpolate (or: interpolation z-value)
#
# Outputs
#   knotsLat: The knots along the latitude axis
#   knotsZ: The knots along the height axis
#   coeffs: The spline coefficients
#   orderLat: The order of the latitude spline
#   orderZ: The order of the height spline
def __determineInterpolationConstants2D__(lat, z, v, degree=None, ax1=None, ax2=None):
    # If necessary reinterpret data
    v = np.asarray(v)
    
    # Construct spline    
    if degree is None:
        degree = 3

    spline = scp_ip.RectBivariateSpline(lat, z, v, kx=degree, ky=degree)
    
    # Draw spline if required
    if (not (ax1 is None)) and (not (ax2 is None)):
        latRefined = np.linspace(lat[0], lat[-1], len(lat) * 100)
        zRefined = np.linspace(z[0], z[-1], len(z) * 100)
        
        # Draw everything with steps in latitude
        for i in range(0, len(lat)):
            vRefined = spline(lat[i], zRefined)
            label = 'lat = ' + str(round(lat[i], 2))
            ax1.plot(zRefined, vRefined[0, :], label=label)
            ax1.plot(z, v[i, :], 'r+')
            
        ax1.legend()
        
        # Draw everything with steps in height
        for i in range(0, len(z)):
            vRefined = spline(latRefined, z[i])
            label = 'z = ' + str(round(z[i], 2))
            ax2.plot(latRefined, vRefined[:, 0], label=label)
            ax2.plot(lat, v[:, i], 'r+')
            
        ax2.legend()
        
    # Return the spline parameters
    knots = spline.get_knots()
    return [knots[0], knots[1], spline.get_coeffs(), degree, degree]

# Atmosphere is the main class that should be used to retrieve atmospheric 
# properties. Once instantiated the functions pressure(), density() and 
# temperature() can be used to retrieve information.
# Instantiation is done by providing a 'definition' argument. This is used to,
# in the future, support multiple atmospheric models such that old code stays
# working as intended. At the moment only the 'preliminary' model is 
# implemented. This model is based on the old VIRA model but is supplemented
# with very basic models to determine wind speed. In the future different
# model definitions will be supported.
# Every function accepts height, latitude and solar longitude as argument. As
# of yet only height is considered in the deep atmosphere (below 33 km) and
# only height and latitude are considered in the upper atmosphere (to 100 km).
# In the future the solar longitude argument will probably be implemented as
# well.
#
# Use of this class is as following (roughly):
#   atm = Atmosphere('preliminary')
#   pressure = atm.pressure(11000, 30, 0, False)
#   minPressure, meanPressure, maxPressure = atm.pressure(11000, 30, 0)

class Atmosphere:
    def __init__(self, definition="preliminary"):
        # To be absolute clear where these values come from: These are 
        # generated from within the __determineInterpolationConstants__ 
        # function using the various data sources in the .csv files
        if definition == 'preliminary':
            self.constants = atm_const.AtmospherePreliminary()
        else:
            raise ValueError("Unknown model identifier")
            
    def __checkAndModifyParameters__(self, height, latitude, solarLongitude, tkc):
        # If any of the variables is an array, convert all the other scalars
        # to a similarly sized array
        arraySize = 0
        
        if util.isArray(height):
            arraySize = len(height)
        elif util.isArray(latitude):
            arraySize = len(latitude)
        elif util.isArray(solarLongitude):
            arraySize = len(solarLongitude)
            
        if arraySize != 0:
            if not util.isArray(height):
                height = np.asarray([height] * arraySize)
            else:
                height = np.asarray(height)
            
            if not util.isArray(latitude):
                latitude = np.asarray([latitude] * arraySize)
            else:
                latitude = np.asarray(latitude)
            
            if not util.isArray(solarLongitude):
                solarLongitude = np.asarray([solarLongitude] * arraySize)
            else:
                solarLongitude = np.asarray(solarLongitude)
                
        latitude = abs(latitude)
        
        if arraySize != 0:
            # Cap latitude values (do not interpolate above or below the
            # known bounds)
            for i in range(0, len(latitude)):
                if latitude[i] < tkc[1][0]:
                    latitude[i] = tkc[1][0]
                elif latitude[i] > tkc[1][-1]:
                    latitude[i] = tkc[1][-1]
                
            for i in range(0, len(height)):
                if height[i] < 0 or height[i] > tkc[0][-1]:
                    raise ValueError("Height is outside of the interpolation bounds")
                    
        else:
            if height < 0 or height > tkc[0][-1]:
                raise ValueError("Height is outside of the interpolation bounds")
                
        return height, latitude, solarLongitude
                
    def __interpolateParameters__(self, height, latitude, tkcDeep, tkcUpper):
        # Preallocate result array and start looping through all values
        isScalar = not util.isArray(height)
        results = []
        
        if isScalar:
            results = [0]
            height = [height]
            latitude = [latitude]
        else:
            results = np.zeros(height.shape)
        
        for i in range(0, len(height)):
            # Check where the height is with respect to the interpolation limits
            if height[i] <= tkcDeep[0][-1]:
                results[i] = scp_ip.splev(height[i], tkcDeep)
            elif height[i] >= tkcUpper[0][0]:
                results[i] = scp_ip.bisplev(height[i], latitude[i], tkcUpper)
            else:
                # Interpolate between the lower and upper interpolating functions (do
                # so linearly for now)
                low = scp_ip.splev(tkcDeep[0][-1], tkcDeep)
                high = scp_ip.bisplev(tkcUpper[0][0], latitude[i], tkcUpper)
                
                results[i] = low + (high - low) * (height[i] - tkcDeep[0][-1]) / \
                    (tkcUpper[0][0] - tkcDeep[0][-1])
                    
        if isScalar:
            return results[0]
            
        return results
        
    def pressure(self, height, latitude, solarLongitude, includeUncertainty=True):
        height, latitude, solarLongitude = \
            self.__checkAndModifyParameters__(height, latitude, solarLongitude, 
                                              self.constants.tkcUpperPres)
        
        mean = self.__interpolateParameters__(height, latitude,
                                              self.constants.tkcDeepPres,
                                              self.constants.tkcUpperPres)
        
        if includeUncertainty == True:
            # Note: uncertainty is given in percentage
            uncertainty = scp_ip.splev(height, self.constants.tkcUncertaintyPres)
            return [mean - uncertainty * mean, mean, mean + uncertainty * mean]
            
        return mean
        
    def density(self, height, latitude, solarLongitude, includeUncertainty=True):
        height, latitude, solarLongitude = \
            self.__checkAndModifyParameters__(height, latitude, solarLongitude, 
                                          self.constants.tkcUpperDens)
        
        mean = self.__interpolateParameters__(height, latitude,
                                              self.constants.tkcDeepDens,
                                              self.constants.tkcUpperDens)
        
        if includeUncertainty == True:
            # Note: uncertainty is given in percentage
            uncertainty = scp_ip.splev(height, self.constants.tkcUncertaintyDens)
            return [mean - (uncertainty * mean), mean, mean + (uncertainty * mean)]
        
        return mean
        
    def temperature(self, height, latitude, solarLongitude, includeUncertainty=True):
        height, latitude, solarLongitude = \
            self.__checkAndModifyParameters__(height, latitude, solarLongitude, 
                                          self.constants.tkcUpperTemp)
        
        mean = self.__interpolateParameters__(height, latitude,
                                              self.constants.tkcDeepTemp,
                                              self.constants.tkcUpperTemp)
        
        if includeUncertainty == True:
            # Note: Uncertainty is given as an absolute value
            uncertainty = scp_ip.splev(height, self.constants.tkcUncertaintyTemp)
            return [mean - uncertainty, mean, mean + uncertainty]
            
        return mean
        
    def velocityMeridional(self, height, latitude, solarLongitude, includeUncertainty=True):
        height, latitude, solarLongitude = \
            self.__checkAndModifyParameters__(height, latitude, solarLongitude, 
                                          self.constants.tkcMeridionalMean)
        
        if includeUncertainty == True:
            return [
                scp_ip.splev(height, self.constants.tkcMeridionalMin),
                scp_ip.splev(height, self.constants.tkcMeridionalMean),
                scp_ip.splev(height, self.constants.tkcMeridionalMax)
            ]
            
        return scp_ip.splev(height, self.constants.tkcMeridionalMean)
        
    def velocityZonal(self, height, latitude, solarLongitude, includeUncertainty=True):
        height, latitude, solarLongitude = \
            self.__checkAndModifyParameters__(height, latitude, solarLongitude, 
                                          self.constants.tkcZonalMean)
        
        if includeUncertainty == True:
            return [
                scp_ip.splev(height, self.constants.tkcZonalMin),
                scp_ip.splev(height, self.constants.tkcZonalMean),
                scp_ip.splev(height, self.constants.tkcZonalMax)            
            ]
            
        return scp_ip.splev(height, self.constants.tkcZonalMean)
        
def __printSplineCoefficients__(variableName, variable, valuesPerLine=5, baseTab=1):
    if not isinstance(variable, list):
        if not hasattr(variable, "__len__"):
            raise ValueError("Expected variable 'variable' to be a list")
    
    tabs = ''
    
    for i in range(0, baseTab):
        tabs += '\t'
    
    msg = tabs + variableName + " = [ "
    for iVariable in range(0, len(variable)):
        if util.isArray(variable[iVariable]):
            numValues = len(variable[iVariable])
            numComplete = int((numValues - (numValues % valuesPerLine)) / valuesPerLine)
            numRemaining = numValues % valuesPerLine
            
            if (numRemaining == 0):
                numRemaining = valuesPerLine
                numComplete -= 1
            
            if iVariable != 0:
                msg += tabs + "[ "
            else:
                msg += "[ "
            
            for iComplete in range(0, numComplete):
                for iSub in range(0, valuesPerLine):
                    msg += str(variable[iVariable][iComplete * valuesPerLine + iSub]) + ", "
                    
                msg += "\n\t" + tabs
            
            for iRemaining in range(0, numRemaining - 1):
                msg += str(variable[iVariable][numComplete * valuesPerLine + iRemaining]) + ", "
                
            msg += str(variable[iVariable][numComplete * valuesPerLine + numRemaining - 1]) + " ]"
        else:
            msg += tabs + str(variable[iVariable])
        
        if iVariable != len(variable) - 1:
            msg += ",\n\t"
    
    msg += " ]"
    
    print(msg)
    
# __determineInterpolationConstants__ is the function determining all
# interpolation coefficients needed to interpolate data within the atmosphere,
# if required the obtained values can be visually inspected.
def __determineInterpolationConstants__(plot=False, valuesPerLine=6):
    # Start by loading the deep atmosphere data. The assumption is that the
    # leftmost header contains altitudes and the topmost header contains, in
    # order: temperature [K], pressure [mbar] and density [kg/m3]
    data = np.loadtxt("./data/atmosphere/VIRADeepAtmosphereProperties.csv",
                      delimiter=';', skiprows=1)
    
    if len(data.shape) != 2 or data.shape[1] != 4:
        raise ValueError('Expected four columns in a 2D deep atmosphere model file')
    
    # retrieve height and make sure it is ascending
    orderedData = np.zeros([data.shape[0] - 1, 4])
    
    if util.isDescending(data[1:, 0]):
        for i in range(0, 4):
            orderedData[:, i] = list(reversed(data[1:, i]))
    else:
        for i in range(0, 4):
            orderedData[:, i] = data[1:, i]
        
    if not util.isAscending(orderedData[:, 0]):
        raise ValueError('Expected height values to be descending or ascending')
        
    # interpolate the various parameters (note that pressure is converted to SI
    # units: from bar to Pa)
    __printSplineCoefficients__('self.tkcDeepTemp', 
        __determineInterpolationConstants1D__(orderedData[:, 0] * 1e3, orderedData[:, 1]))
    __printSplineCoefficients__('self.tkcDeepPres',
        __determineInterpolationConstants1D__(orderedData[:, 0] * 1e3, orderedData[:, 2] * 1e5))
    __printSplineCoefficients__('self.tkcDeepDens',
        __determineInterpolationConstants1D__(orderedData[:, 0] * 1e3, orderedData[:, 3]))
    
    # Now load the upper atmosphere data. The assumption is now that the 
    # leftmost header contains altitude data. The rest of the data comes in
    # grouped columns of 3. Each group is for a given latitude, described at
    # the top. and contains temperature, pressure and density data 
    # consecutively
    data = np.genfromtxt("./data/atmosphere/VIRAUpperAtmosphereProperties.csv",
                      delimiter=';')
    
    if len(data.shape) != 2 or data.shape[1] < 4:
        raise ValueError("Expected a 2D upper atmosphere model file with at least one contained dataset")
        
    if (data.shape[1] - 1) % 3 != 0:
        raise ValueError("Expected latitude data to come in sets of 3")
    
    # Preallocate arrays to store data in    
    numHeightValues = data.shape[0] - 2
    numLatitudeValues = int((data.shape[1] - 1) / 3)
    
    height = np.zeros(numHeightValues)
    latitude = np.zeros(numLatitudeValues)
    
    temp = np.zeros([numHeightValues, numLatitudeValues])
    pressure = np.zeros([numHeightValues, numLatitudeValues])
    density = np.zeros([numHeightValues, numLatitudeValues])
    
    # Extract data
    for iLatitude in range(0, numLatitudeValues):
        latitude[iLatitude] = float(data[0, 1 + iLatitude * 3])
        
    for iHeight in range(0, numHeightValues):
        height[iHeight] = float(data[2 + iHeight, 0])
        
        for iLatitude in range(0, numLatitudeValues):
            temp[iHeight, iLatitude] = float(data[2 + iHeight, 1 + iLatitude * 3])
            pressure[iHeight, iLatitude] = float(data[2 + iHeight, 2 + iLatitude * 3])
            density[iHeight, iLatitude] = float(data[2 + iHeight, 3 + iLatitude * 3])
            
    # Check if the data is ordered correctly
    if util.isDescending(height):
        height = np.flipud(height)
        temp = np.flipud(temp)
        pressure = np.flipud(pressure)
        density = np.flipud(density)
        
    if not util.isAscending(height):
        raise ValueError("Upper atmosphere height is neither ascending nor descending")

    __printSplineCoefficients__('self.tkcUpperTemp',
        __determineInterpolationConstants2D__(height * 1e3, latitude, temp))
    __printSplineCoefficients__('self.tkcUpperPres',
        __determineInterpolationConstants2D__(height * 1e3, latitude, pressure * 1.0e5))
    __printSplineCoefficients__('self.tkcUpperDens',
        __determineInterpolationConstants2D__(height * 1e3, latitude, density))

    # Load the wind speed data
    data = np.genfromtxt("./data/atmosphere/VIRAWindSpeeds.csv",
                          dtype=str, delimiter=';')
    
    if len(data.shape) != 2 or data.shape[1] != 5 or data.shape[0] < 3:
        raise ValueError("Expected a 2D wind speed file with exactly 5 columns and at least 3 columns")
    
    # Look for the zonal and meridional keywords and determine number of rows
    iZonal = 0
    iMeridional = 0
    
    for i in range(1, data.shape[0]):
        if isinstance(data[i][0], str):
            if data[i][0].lower() == 'zonal':
                iZonal = i
                if iMeridional != 0:
                    break
            elif data[i][0].lower() == 'meridional':
                iMeridional = i
                if iZonal != 0:
                    break
                
    if iZonal == 0 or iMeridional == 0:
        raise VaueError("Failed to find the 'zonal' and/or 'meridional' keywords")
            
    # determine iteration indices
    numZonal = 0
    numMeridional = 0
    endZonal = 0
    endMeridional = 0
    
    if iZonal > iMeridional:
        numZonal = data.shape[0] - iZonal
        endZonal = data.shape[0]
        numMeridional = iZonal - iMeridional
        endMeridional = iZonal
    else:
        numMeridional = data.shape[0] - iMeridional
        endMeridional = data.shape[0]
        numZonal = iMeridional - iZonal
        endZonal = iMeridional
        
    # preallocate arrays
    heightZonal = np.zeros(numZonal)
    minZonal = np.zeros(numZonal)
    meanZonal = np.zeros(numZonal)
    maxZonal = np.zeros(numZonal)
    
    heightMeridional = np.zeros(numMeridional)
    minMeridional = np.zeros(numMeridional)
    meanMeridional = np.zeros(numMeridional)
    maxMeridional = np.zeros(numMeridional)
    
    # retrieve information from dataset
    for i in range(iZonal, endZonal):
        iAssign = i - iZonal
        heightZonal[iAssign] = float(data[i][1])
        minZonal[iAssign] = float(data[i][2])
        maxZonal[iAssign] = float(data[i][3])
        meanZonal[iAssign] = float(minZonal[iAssign] + maxZonal[iAssign]) / 2.0
        
    for i in range(iMeridional, endMeridional):
        iAssign = i - iMeridional
        heightMeridional[iAssign] = float(data[i][1])
        minMeridional[iAssign] = float(data[i][2])
        maxMeridional[iAssign] = float(data[i][3])
        meanMeridional[iAssign] = float(minMeridional[iAssign] + maxMeridional[iAssign]) / 2.0
        
    # Check if the data is ordered correctly
    if util.isDescending(heightZonal):
        heightZonal = np.flipud(heightZonal)
        minZonal = np.flipud(heightZonal)
        meanZonal = np.flipud(meanZonal)
        maxZonal = np.flipud(maxZonal)
        
    if util.isDescending(heightMeridional):
        heightMeridional = np.flipud(heightMeridional)
        minMeridional = np.flipud(minMeridional)
        meanMeridional = np.flipud(meanMeridional)
        maxMeridional = np.flipud(maxMeridional)
        
    if not util.isAscending(heightZonal):
        raise ValueError("Zonal wind speed height should be purely ascending or descending")
        
    if not util.isAscending(heightMeridional):
        raise ValueError("Meridional wind speed height should be purely ascending or descending")
        
    __printSplineCoefficients__("self.tkcZonalMin", 
        __determineInterpolationConstants1D__(heightZonal * 1e3, minZonal, degree=2))
    __printSplineCoefficients__("self.tkcZonalMean",
        __determineInterpolationConstants1D__(heightZonal * 1e3, meanZonal, degree=2))
    __printSplineCoefficients__("self.tkcZonalMax",
        __determineInterpolationConstants1D__(heightZonal * 1e3, maxZonal, degree=2))
    
    __printSplineCoefficients__("self.tkcMeridionalMin",
        __determineInterpolationConstants1D__(heightMeridional * 1e3, minMeridional, degree=2))
    __printSplineCoefficients__("self.tkcMeridionalMean",
        __determineInterpolationConstants1D__(heightMeridional * 1e3, meanMeridional, degree=2))
    __printSplineCoefficients__("self.tkcMeridionalMax",
        __determineInterpolationConstants1D__(heightMeridional * 1e3, maxMeridional, degree=2))
    
    # Open and process uncertainty files
    files = ['./data/atmosphere/VIRAUncertaintyTemperature.csv',
             './data/atmosphere/VIRAUncertaintyPressure.csv',
             './data/atmosphere/VIRAUncertaintyDensity.csv']
    variables = ['self.tkcUncertaintyTemp',
                 'self.tkcUncertaintyPres',
                 'self.tkcUncertaintyDens']
    divFactor = [1, 100, 100]
                 
    for iVar in range(0, len(variables)):
        data = np.genfromtxt(files[iVar], delimiter=';', skip_header=1)
        
        __printSplineCoefficients__(variables[iVar],
            __determineInterpolationConstants1D__(data[:,0] * 1e3, data[:,1] / divFactor[iVar], degree=1))
    
# __testPressure__ is a simple functions plotting the pressure variation
# throughout the lower atmosphere and throughout the upper atmopshere for 
# manual inspection
def __testPressure__():
    _, ax = plt.subplots(1, 1)
    
    atm = Atmosphere()
    zDeep = np.linspace(atm.constants.tkcDeepPres[0][0], atm.constants.tkcDeepPres[0][-1], 1000)
    presDeep = atm.pressure(zDeep, [0] * len(zDeep), 0)
    ax.plot(presDeep[0] / 1e5, zDeep / 1e3, 'g')
    ax.plot(presDeep[1] / 1e5, zDeep / 1e3, 'b')
    ax.plot(presDeep[2] / 1e5, zDeep / 1e3, 'r')
    ax.legend(['min', 'mean', 'max'])
    ax.set_xlabel('Pressure [bar]')
    ax.set_ylabel('Height [km]')
    
def __testTemperature__():
    _, ax = plt.subplots(1, 1)
    
    atm = Atmosphere()
    zDeep = np.linspace(atm.constants.tkcDeepTemp[0][0], atm.constants.tkcUpperTemp[0][-1], 1000)
    tempDeep = atm.temperature(zDeep, [0] * len(zDeep), 0)
    ax.plot(tempDeep[0], zDeep / 1e3, 'g')
    ax.plot(tempDeep[1], zDeep / 1e3, 'b')
    ax.plot(tempDeep[2], zDeep / 1e3, 'r')
    ax.legend(['min', 'mean', 'max'])
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Height [km]')
    
def __testDensity__():
    _, ax = plt.subplots(1, 1)
    
    atm = Atmosphere()
    zDeep = np.linspace(atm.constants.tkcDeepDens[0][0], atm.constants.tkcDeepDens[0][-1], 1000)
    densDeep = atm.density(zDeep, [0] * len(zDeep), 0)
    ax.plot(densDeep[0], zDeep / 1e3, 'g')
    ax.plot(densDeep[1], zDeep / 1e3, 'b')
    ax.plot(densDeep[2], zDeep / 1e3, 'r')
    ax.legend(['min', 'mean', 'max'])
    ax.set_xlabel('Density [kg/m3]')
    ax.set_ylabel('Height [km]')
    
def __testVelocity__():
    _, [ax1, ax2] = plt.subplots(1, 2)
    
    atm = Atmosphere()
    zDeep = np.linspace(atm.constants.tkcZonalMean[0][0], atm.constants.tkcZonalMean[0][-1], 1000)
    zonal = atm.velocityZonal(zDeep, 0, 0)
    print('zonal =', np.asarray(zonal).shape)
    ax1.plot(zonal[0], zDeep / 1e3, 'g')
    ax1.plot(zonal[1], zDeep / 1e3, 'b')
    ax1.plot(zonal[2], zDeep / 1e3, 'r')
    ax1.legend(['min', 'mean', 'max'], loc='best')
    ax1.set_xlabel('Zonal velocity [m/s]')
    ax1.set_ylabel('Height [km]')
    
    meridional = atm.velocityMeridional(zDeep, 0, 0)
    print('meridional =', np.asarray(zonal).shape)
    ax2.plot(meridional[0], zDeep / 1e3, 'g')
    ax2.plot(meridional[1], zDeep / 1e3, 'b')
    ax2.plot(meridional[2], zDeep / 1e3, 'r')
    ax2.legend(['min', 'mean', 'max'], loc='best')
    ax2.set_xlabel('Meridional velocity [m/s]')
    ax2.set_ylabel('Height [km]')
    
def __testScalar__():
    atm = Atmosphere()
    atm.pressure(0, 5, 5)

#__testScalar__() 
#__testPressure__()
#__testDensity__()
#__testTemperature__()   
#__testVelocity__()
#__determineInterpolationConstants__()