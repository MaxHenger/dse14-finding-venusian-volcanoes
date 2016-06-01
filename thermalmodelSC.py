# -*- coding: utf-8 -*-
"""
Created on Mon May 30 14:48:46 2016

@author: Chaggai
"""
import numpy as np 
import SCpower 
import HeatTransferFunctions as trans
import matplotlib.pyplot as plt

def Spacecraftthermal(distsun, abs_solar_SC, abs_IR_SC, A_SC_sol, A_SC_plan, A_total, albedo_venus, R_venus, R_orbit, Tmin, torbit, teclipse, heatcap, mass_SC, P_ins, eff_ins,heater_P,distsun_e):
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
    Q_abs_IR =  flux_IR*abs_IR_SC*A_SC_plan
    
    Qabs_tot =Q_abs_IR+Q_abs_sol+ Q_abs_alb+P_ins*(1.-eff_ins)
    
    Q_emit = Qabs_tot
    T_SC =  (Q_emit/(abs_IR_SC*5.670373*10**-8*A_total))**(1./4.)
    
    Q_emit_min = abs_IR_SC*5.670373*10**-8*A_total*Tmin**4
    Q_heat = Q_emit_min-Q_abs_IR-+P_ins*(1.-eff_ins)
    
    #interplanetary 
    Q_abs_solinter = solfluxinter*abs_solar_SC*A_SC_sol+320.
    T_SC_inter =  (Q_abs_solinter/(abs_IR_SC*5.670373*10**-8*A_total))**(1./4.)-273.15
    print T_SC_inter
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
            Q_out = abs_IR_SC*5.670373*10**-8*A_total*T**4
            
            Qbal = Qin-Q_out
            dT = Qbal*1/(mass_SC*heatcap)
            Tnew = T+dT
            tempplot.append(Tnew-273.15)
    plt.plot(timeplot,tempplot)
    
    return(T_SC,Q_heat)
    
if __name__=="__main__":
    abs_IR_SC = 0.21#0.3*0.8 #http://www.sheldahl.com/documents/RedBook.pdf  aluminium+1080 layer
    abs_solar_SC = 0.205 #0.22*.85
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
    heatcap = 600.
    mass_SC = 1500.
    P_ins = 500.
    P_eff = 0.2
    heater_P = 000.
    de = 149597870700
    print Spacecraftthermal(distsun,abs_solar_SC, abs_IR_SC, A_SC_sol, A_SC_plan, A_total, albedo_venus, R_venus, R_orbit,Tmin, torbit, teclipse, heatcap, mass_SC, P_ins, P_eff,heater_P,de)