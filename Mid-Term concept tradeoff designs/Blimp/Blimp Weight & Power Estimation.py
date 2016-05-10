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
Pwr_ther = 700 #from a fridge

Mpl = 100.


def blimpweight(Mblp, x):
    Q = Mblp/rho_at
    r = sp.special.cbrt(3*Q/4./np.pi)
    A = 4.*np.pi*r**2
    if x == 'solar':
        Asolar = 0.86*A
        Msolar = Asolar*rho_solar
        Mbat = rho_bat*t_ecl
    elif x == 'rtg':
        Asolar = np.zeros(np.size(Mblp))
        Msolar = 230. #for 1200 W RTG
        Mbat = 0
    t = Pat*r/2./sig_kav
    Mskin = rho_kavler*t*A*sf
    Mh2 = rho_h2*Q

    Mttc = ratio_ttc* Mblp
    Madcs = ratio_adcs * Mblp
    Mstruc = ratio_struc * Mblp

    Mproptherm = Mblp - Mskin - Mh2 - Msolar - Mttc - Madcs - Mstruc - Mbat - Mpl
    return Mproptherm, Mskin, Mh2, Msolar, Mttc, Madcs, Mstruc, Mbat, r, Q, Asolar

def blimppwr(Mblp, x):
    if x == 'solar':
        Pwrblp = pwrrho*blimpweight(Mblp, x)[-1]
        Pwr_pl_prop_therm = Pwrblp/2 - pwrratio_pwr_ttc_adcs* Pwrblp
        Pwr_pl_prop = np.zeros(np.size(Mblp))
    elif x == 'rtg':
        Pwrblp = np.ones(np.size(Mblp))*1200.
        Pwr_pl_prop_therm = np.zeros(np.size(Mblp))
        Pwr_pl_prop = Pwrblp - pwrratio_pwr_ttc_adcs* Pwrblp - Pwr_ther
    return Pwr_pl_prop_therm, Pwrblp, Pwr_pl_prop

x = 'rtg'

Mblp = np.arange(500, 2010, 10)

Mpropthermlist = blimpweight(Mblp, x)[0]
Msolar = blimpweight(Mblp, x)[3]
r = blimpweight(Mblp, x)[-3]
Q = blimpweight(Mblp, x)[-2]
Asolar = blimpweight(Mblp, x)[-1]

Pwr_pl_prop_thermlist = blimppwr(Mblp, x)[0]
Pwrblp = blimppwr(Mblp, x)[1]


print blimppwr(1500, x)[2]

plt.plot(Mblp, Pwr_pl_prop_thermlist, 'r')
plt.title('Power payload, prop, thermal vs Blimp weight')
plt.show()
plt.plot(Mblp, Pwrblp)
plt.title('Blimp power vs weight')
plt.show()
plt.plot(Mblp, Mpropthermlist, 'r')
#plt.ylim([0, 250])
plt.title('Prop, thermal Weight vs Blimp Weight')
plt.show()

#plt.plot(Mblp, Msolar)
#plt.title('Sloar weight vs Blimp weight')
#plt.show()
#plt.plot(Asolar, Q)
#plt.title('Volume vs Area')
#plt.show()






