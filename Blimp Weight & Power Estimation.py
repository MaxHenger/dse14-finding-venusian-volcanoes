# -*- coding: utf-8 -*-
"""
Created on Mon May 09 10:57:17 2016

@author: Yuyang
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

rho_at = 3.
g = 8.8
sf = 1.5
rho_kavler = 1440. #kg/m3
sig_kav = 3600000000. #Pa
Pat = 8. * 100000 #Pa
rho_h2 = 0.14 #kg/m3
rho_solar = 0.35 #kg/m2

ratio_ttc = 0.026 #weight percentage of blimp mass
ratio_adcs = 0.034 #weight percentage of blimp mass
ratio_struc = 0.282 #weight percentage of blimp mass
rho_bat = 3.77 #kg/h
t_ecl = 13.46*sf

pwrrho = 2.26 #W/m2
pwrratio_pwr_ttc_adcs = 0.17 #Power ratio of power, TTC, ADCS

Mpl = 50.

def blimpweight(Mblp):
    Q = Mblp/rho_at
    r = sp.special.cbrt(3*Q/4./np.pi)
    A = 4.*np.pi*r**2
    Asolar = 0.86*A
    t = Pat*r/2./sig_kav
    Mskin = rho_kavler*t*A*sf
    Mh2 = rho_h2*Q
    Msolar = Asolar*rho_solar
    Mttc = ratio_ttc* Mblp
    Madcs = ratio_adcs * Mblp
    Mstruc = ratio_struc * Mblp
    Mbat = rho_bat*t_ecl
    Mproptherm = Mblp - Mskin - Mh2 - Msolar - Mttc - Madcs - Mstruc - Mbat - Mpl
    return Mproptherm, Mskin, Mh2, Msolar, Mttc, Madcs, Mstruc, Mbat, r, Q, Asolar

def blimppwr(Mblp):
    Pwrblp = pwrrho*blimpweight(Mblp)[-1]
    Pwr_pl_prop_therm = Pwrblp/2 - pwrratio_pwr_ttc_adcs* Pwrblp
    return Pwr_pl_prop_therm, Pwrblp

Mblplist = np.arange(800, 1210, 10)

Mpropthermlist = blimpweight(Mblplist)[0]
Msolar = blimpweight(Mblplist)[3]
r = blimpweight(Mblplist)[-3]
Q = blimpweight(Mblplist)[-2]

Pwr_pl_prop_thermlist = blimppwr(Mblplist)[0]
Pwrblp = blimppwr(Mblplist)[1]



a = blimpweight(1000)

plt.plot(Mblplist, Pwr_pl_prop_thermlist, 'r')
plt.plot(Mblplist, Pwrblp)

#plt.plot(Mblplist, Mpropthermlist, 'r')
#plt.ylim([0, 250])
#
#plt.plot(Mblplist, Msolar)
#
#plt.show()
#plt.plot(Q, r)

plt.show()






