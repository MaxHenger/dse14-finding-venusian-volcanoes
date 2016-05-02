"""
Created on Mon May 02 14:22:24 2016

@author: Mathijs

Solar power calculator
"""

from scipy import interpolate
import scipy as sp


# Insert altitude in km, panel area in meters squared and angle of incidence in radians
# NOT TO BE USED ABOVE 60 KM!!
def SolarPower(h,A_Panel,incidence):
    
    # Define altitude, solar irradiance and efficiency from literature study for altitude 0-60 km (http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20150016298.pdf)
    altlist = sp.arange(0,65,5)
    Irradiance = [8.7,22.0,37.4,45.6,58.5,71.9,85.1,99.0,112.2,181.4,256.0,404.4,559.2]
    Efficiency = [0.33,0.84,1.43,1.75,2.24,2.75,3.26,3.79,4.29,6.94,9.79,15.47,21.39]   
    
    # Interpolate irradiance and efficiency for altitude provided
    Irr_inter = sp.interpolate.InterpolatedUnivariateSpline(altlist,Irradiance)
    Eff_inter = sp.interpolate.InterpolatedUnivariateSpline(altlist,Efficiency)
    Irr_inter = sp.interpolate.InterpolatedUnivariateSpline(altlist,Irradiance)

#    print Irr_inter(h)
#    print sp.sin(incidence)
#    print Eff_inter(h)

    # Determine solar power generated in W
    P_solar = Irr_inter(h) * (Eff_inter(h)/100.) * sp.sin(incidence)
    
    return P_solar

P_solar = SolarPower(60,1,0.5*sp.pi)
print P_solar