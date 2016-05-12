# -*- coding: utf-8 -*-
"""
Created on Tue May 10 13:40:28 2016

@author: Jeffrey
"""

import math

class massSettings:
    def __init__(self):
        self.rho = 0 # set this
        self.h = 0 # set this
        self.t = 0.004
        self.t_sheet = 0.002
        self.b = 0 # set this
        self.c = 0 # set this
        self.A_stiff = 0.0004
        self.x_stiff = 0.1
        self.x_ribt = 3.0
        self.R = 0 # set this
        self.t_f = 0.002
        self.L = 0 # set this

##spar mass
def sparmass(settings):
    m_spar = 2*settings.rho*settings.h*settings.t*settings.b
    return m_spar

##sheet mass
def sheetmass(settings):
    m_sheet = 2*settings.rho*settings.t_sheet*settings.c*settings.b
    return m_sheet

##stiffener mass
def stiffmass(settings):
    m_stiff = settings.c/settings.x_stiff * settings.rho*settings.A_stiff * settings.b
    return m_stiff

##rib mass
def ribmass(settings):
    m_rib = settings.b/settings.x_rib * settings.rho*settings.h*settings.t*settings.c
    return m_rib

##wing mass
def wingmass(settings):
    m_spar = sparmass()
    m_sheet = sheetmass()
    m_stiff = stiffmass()
    m_rib = ribmass()
    m_w = m_spar + m_sheet + m_stiff + m_rib 
    return m_w

##tail mass
def tailmass(settings):
    m_w = wingmass()
    m_t = m_w/4
    return m_t

####fuselage thickness
##def fusethick():
##    t_f = p*R/sigma
##    return t_f

##fuselage mass
def fusemass(settings):
    m_f = 2*math.pi*settings.R*settings.t_f*settings.L*settings.rho
    return m_f

##propulsion system mass
def propmass(settings):
    m_w = wingmass()
    m_p = m_w/4.0
    return m_p

def allmass(settings):
    return wingmass(settings) + tailmass(settings) + \
        fusemass(settings) + propmass(settings)