# -*- coding: utf-8 -*-
"""
Created on Wed May 25 10:30:50 2016

@author: Chaggai
"""
import numpy as np

def solarflux(r):
    Lsun=382.8*10**24 #Watts
    Solfluxloc = Lsun/(4.*np.pi*r**2)
    return(Solfluxloc)

def SCsolarr(Puseday, Pusenight, meanicidenceangle,rorbit, torbit, teclipse, eff_PV, eff_battery, eff_EPS):
    tday = torbit-teclipse
    effday = eff_EPS
    effeclip = eff_battery*eff_EPS
    Eneccday = Puseday/effday*(tday/3600.) #Wh
    Enecceclip = Pusenight/effeclip*(teclipse/3600.)#wh
    Etot =Eneccday+Enecceclip #in Watt hour
    PprovidePV = Etot/(tday/3600.) #watts to be provided during day
    solflux = solarflux(rorbit) #rorbit is distance from sun in meters
    Area=PprovidePV/(solflux*np.cos(np.deg2rad(meanicidenceangle)))/eff_PV
    return(Area)    
    
    
def batterysize(Pusenight,teclipse, DOD, diseff, specenergy,energdens):
    Ebat = Pusenight*(teclipse/3600.)/(DOD*diseff)
    Mass_bat = Ebat/specenergy
    vol_bat = Ebat/energdens
    return(Ebat, Mass_bat, vol_bat)
    
if __name__=="__main__":
    smaearth = 149.6*10**6*1000. #m
    smavenus = 108.21*10**6*1000. #m
    eff_PV = 0.29 #http://www.spectrolab.com/DataSheets/cells/2015%20XTJ%20CIC%20Datsheet.pdf
    degradeyear = 0.98
    missionlife = 5 #years
    eff_PV = eff_PV*degradeyear**missionlife
    eff_battery = 0.6 
    eff_EPS = 0.8 #ADSEE 1 reader
    Puseday = 2000.#Watts
    Pusenight = 1500 #watts
    Porbit = 12000.#Seconds
    teclipse = 0.4*Porbit#seconds
    meanicidenceangle = 30.#degrees mean so average angle during day part of one orbit
    
    
    #battery specifications
    DOD = 0.25 #http://www.saftbatteries.com/force_download/li_ion_battery_life__TechnicalSheet_en_0514_Protected.pdf
    diseff = 0.9 
    specenergy = 170. #wh/kg taken from saft lion datasheet
    energdens = 250. #Wh/l
    print SCsolarr(Puseday,Pusenight, meanicidenceangle, smavenus, Porbit, teclipse, eff_PV, eff_battery, eff_EPS)
    print batterysize(Pusenight,teclipse, DOD, diseff, specenergy,energdens)
    
    
    #size mass cables when lines are known