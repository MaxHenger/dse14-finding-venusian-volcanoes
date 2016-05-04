# -*- coding: utf-8 -*-
"""
Created on Tue May  3 12:29:47 2016

@author: MaxHenger
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as scp_ip

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
        spline = scp_ip.splrep(z, v)
    else:
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
    spline = []
    
    if degree is None:
        spline = scp_ip.RectBivariateSpline(lat, z, v)
    else:
        spline = scp_ip.RectBivariateSpline(lat, z, v, kx=degree, ky=degree)
    
    # Draw spline if required
    if not (ax is None):
        latRefined = np.linspace(lat[0], lat[-1], len(lat) * 100)
        zRefined = np.linspace(z[0], z[-1], len(z) * 100)
        
        # Draw everything with steps in latitude asd asd
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

class Atmosphere:
    def __init__(self):
        self.values = 0
        
# __determineInterpolationConstants__ is the function determining all
# interpolation coefficients needed to interpolate data within the atmosphere,
# if required the obtained values can be visually inspected.
def __determineInterpolationConstants__(file, plot=False):
    