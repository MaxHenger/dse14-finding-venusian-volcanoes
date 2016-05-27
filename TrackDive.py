# -*- coding: utf-8 -*-
"""
Created on Fri May 27 11:58:10 2016

@author: MaxHenger
"""

import Atmosphere
import matplotlib.pyplot as plt

# Analysis_TrackDive is a function that will perform some simple plotting
# for a preliminary investigation into the effect of changing parameters
# such as the required acceleration, changing the top and bottom height
# and other parameters
def Analysis_TrackDive(lt, lb):
    atm = Atmosphere.Atmosphere("preliminary")
    