# -*- coding: utf-8 -*-
"""
Created on Mon May 30 14:48:46 2016

@author: Chaggai
"""
import numpy as np 
import SCpower 
import HeatTransferFunctions as trans
import matplotlib.pyplot as plt

def Spacecraftthermal(distsun, abs_solar_SC, e_IR_SC, A_SC_sol, A_SC_plan, A_total, albedo_venus, R_venus, R_orbit, torbit, teclipse, heatcap, mass_SC, P_ins, eff_ins,heater_P,distsun_e,A_rad, e_rad):
    #solar fluxincoming
    solflux = SCpower.solarflux(distsun)
    solfluxinter = SCpower.solarflux(distsun_e)
    Q_abs_sol = solflux*abs_solar_SC*A_SC_sol
    
    #albedo from planet incoming flux
    visibility_fac = (R_venus/R_orbit)**2
    flux_alb = solflux*albedo_venus*visibility_fac
    Q_abs_alb = flux_alb*abs_solar_SC*A_SC_plan
    
    #IR from planet temp
    T_backbody_venus = 226.6 #http://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html
    flux_IR = 5.670373*10**-8*T_backbody_venus**4
    Q_abs_IR =  flux_IR*e_IR_SC*A_SC_plan
    
    

    #interplanetary 
    Q_abs_solinter = solfluxinter*abs_solar_SC*A_SC_sol+500.
    T_SC_inter =  (Q_abs_solinter/(5.670373*10**-8*(e_IR_SC*(A_total-A_rad)+e_rad*A_rad)))**(1./4.)-273.15
    #time based simulator
    time = np.arange(torbit)
    timeplot = np.arange(torbit*3)
    tempplot = []
    Tnew = 290.
    tday = torbit-teclipse
    for orbit in range(3):
        for t in time:
            T= Tnew
            if t <=tday:
                Qin = Q_abs_IR+Q_abs_sol+ Q_abs_alb+P_ins*(1.-eff_ins)
            if t>tday:
                Qin = Q_abs_IR+P_ins*(1.-eff_ins)+heater_P
            Q_out = e_IR_SC*5.670373*10**-8*(A_total-A_rad)*T**4+e_rad*5.670373*10**-8*A_rad*T**4
            
            Qbal = Qin-Q_out
            dT = Qbal*1/(mass_SC*heatcap)
            Tnew = T+dT
            tempplot.append(Tnew-273.15)
    plt.plot(timeplot,tempplot)
    Tmax_orb = np.max(tempplot)
    Tmin_orb = np.min(tempplot)    
    return(Tmax_orb, Tmin_orb, T_SC_inter)
    
    
    
if __name__=="__main__":
    distsun = 108.21*10**6*1000.
    A_SC_sol = 3.25*8*1.1
    A_SC_plan =A_SC_sol
    A_total = 4.*3.25*8+2.*3.25**2
    albedo_venus=0.65 #nasa fact sheet
    alt_SC = 112022.5*1000.
    R_venus = 6051.8*1000.
    R_orbit = alt_SC+R_venus
    Tmin=273.15
    torbit = P_SC = 2894.9*60.
    teclipse = 0.46*3600.
    P_eff = 0.2
    heater_P = 000.
    de = 149597870700
    
    
    heatcap = 600.
    mass_SC = 1500.
    P_ins = 500.

    e_IR_SC = 0.12644#2 blankets kapton
    abs_solar_SC = 0.14
    e_rad = 0.85
    A_rad =2.8
    


    Tmax_orb, Tmin_orb, T_SC_inter= Spacecraftthermal(distsun, abs_solar_SC, e_IR_SC, A_SC_sol, A_SC_plan, A_total, albedo_venus, R_venus, R_orbit, torbit, teclipse, heatcap, mass_SC, P_ins, P_eff,heater_P,de,A_rad, e_rad)
    print Tmax_orb, Tmin_orb, T_SC_inter
        