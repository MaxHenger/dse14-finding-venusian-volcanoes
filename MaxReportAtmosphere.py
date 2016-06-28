#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 00:49:06 2016

@author: MaxHenger
"""

import Atmosphere


import matplotlib.pyplot as plt
import numpy as np

def PlotAtmosphericProperties():
    atm = Atmosphere.Atmosphere()
    height = np.linspace(0, 85e3, 500)
    density = atm.density(height, 0, 0)
    temperature = atm.temperature(height, 1, 0)
    pressure = atm.pressure(height, 0, 0)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(temperature[0], height / 1e3, 'k--')
    ax.plot(temperature[1], height / 1e3, 'k')
    ax.plot(temperature[2], height / 1e3, 'k--')
    
    ax.set_xlabel(r'$T \; [K]$', fontsize=14)
    ax.set_ylabel(r'$h \; [km]$', fontsize=14)
    ax.grid(True)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(np.log10(pressure[0]), height / 1e3, 'k--')
    ax.plot(np.log10(pressure[1]), height / 1e3, 'k')
    ax.plot(np.log10(pressure[2]), height / 1e3, 'k--')
    
    ax.set_xlabel(r'$\mathrm{log}_{10}( p ) \; [Pa]$', fontsize=14)
    ax.set_ylabel(r'$h \; [km]$', fontsize=14)
    ax.grid(True)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(np.log10(density[0]), height / 1e3, 'k--')
    ax.plot(np.log10(density[1]), height / 1e3, 'k')
    ax.plot(np.log10(density[2]), height / 1e3, 'k--')
    
    ax.set_xlabel(r'$\mathrm{log}_{10}( \rho ) \; [kg/m^3]$', fontsize=14)
    ax.set_ylabel(r'$h \; [km]$', fontsize=14)
    ax.grid(True)
    
def PlotAtmosphericVelocities():
    atm = Atmosphere.Atmosphere()
    height = np.linspace(0, 80e3, 500)
    vZonal = atm.velocityZonal(height, 0, 0)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(vZonal[0], height / 1e3, 'k--')
    ax.plot(vZonal[1], height / 1e3, 'k')
    ax.plot(vZonal[2], height / 1e3, 'k--')
    
    ax.set_xlabel(r'$V_\mathrm{zonal} \; [m/s]$', fontsize=14)
    ax.set_ylabel(r'$h \; [km]$', fontsize=14)
    ax.grid(True)
    
#PlotAtmosphericProperties()
#PlotAtmosphericVelocities()