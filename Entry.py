# -*- coding: utf-8 -*-
"""
Created on Tue May 10 18:46:59 2016

This file contains the three entry methods( balistic, gliding and skipping) in
form of three classes. Each class has seperate methods for retrival of various
unique and general information. All contain a method for expected landing as
well as error bars indicating possible landing ranges. 

Ballistic class:

Gliding class:

Skipping class:

In general all functions accept and return values in SI units. An exception is
angles, these are entered in degrees. Also attempt to call the functions with
arrays as much as possible to reduce unneccesary calculations.

The general use of this file is as following:
    
    import Entry
    import Orbit

    orb = Orbit.Orbit('preliminary')
    orb.create(SemiMajor,Eccentricity,Inclination,AscentionAngle)
    orb.place(Location,Time)
    
    bal = Entry.Ballistic(orb)

@author: Julius
"""

# Given a certain orbit, entry method, landing site: calculate time to release
# orbit and release location, decent profile, decent acceleration, decent heat 
# flux, decent heat

#
#
#

import matplotlib.pyplot as plt
import numpy as np
import Utility as util
import Atmosphere

atm = Atmosphere.Atmosphere("preliminary")


class Ballistic:
    """Ballistic entry profile: takes orbit as inital argument """
    def __init__(self,Orbit)
        self.orbit=Orbit
        





