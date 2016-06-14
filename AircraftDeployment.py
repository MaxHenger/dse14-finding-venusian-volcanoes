"""
Created on Sun May 29 18:02:17 2016

@author: Mathijs
"""

import Atmosphere
import numpy as np
import matplotlib.pyplot as plt

"""
Inputs:
alt = initial altitude at which parachute is fully deployed [m]
mass = mass of entire system [kg]
Vel = inital velocity at which the system works [m/s]
C_D_chute = drag coefficient of parachute [-]
A_chute = area of parachute [m^2]
g_venus = Venusian gravitational acceleration [m/2^2]
dt = time step for calculation [s]
mass_shield = mass of heat shield [kg]
t_shield = time at which heat shield is ejected with respect to deployment [s]
lat = lattitude at which lander is deployed [deg]
sol = solar latitude at which lander is deployed [deg]
SF_drag = factor to account for drag created by package [-]a
alt_end = altitude at which the model has to stop [m]

Outputs:
altarray = array with altitudes up to altitude where aircraft is deployed [m]
timearray = array of time [s]
velarray = array with velocity of system [m/s]

"""

def ParachuteTrajectory(alt,mass,Vel,C_D_chute,A_chute,g_venus,dt,mass_shield,t_shield,lat,sol,SF_drag,alt_end):
    
    # Set up atmosphere for the temperature
    atm = Atmosphere.Atmosphere('preliminary')
    
    lat = 0
    sol = 0

    altarray = np.ones(1)*alt    
    timearray = np.zeros(1)
    velarray = np.ones(1)*Vel
    massarray = np.ones(1)*mass
    accarray = np.zeros(1)
    t = 0    
    first = True
    W = mass*g_venus
    
#    for i in range(1):
    while alt > alt_end:
        # Determine average density at altitude, and position [kg/m^3]  
        Density = (atm.density(alt,lat,sol, includeUncertainty=True)[1])
        
        # Remove weight of shield once t_shield is reached, and only do so once
        if t > t_shield and first:
            mass = mass-mass_shield
            first = False

        # Determine drag force
        DragForce = 0.5*C_D_chute*Density*A_chute*Vel**2*SF_drag
        # Determine resultant force
        ResultForce = W-DragForce
        # Calculate acceleration
        acceleration = ResultForce/mass        
        
        # Determine distance covered in time dt
        Dist = Vel*dt + 0.5*acceleration*dt**2
        # Determine velocity for next phase
        Vel = Vel + acceleration*dt
        
        # Determine altitude
        alt = alt - Dist
        # Determine time passed [s]
        t = t + dt

        # Add values to arrays for plotting and output
        altarray = np.append(altarray,alt)
        timearray = np.append(timearray,t)
        velarray = np.append(velarray,Vel)
        massarray = np.append(massarray,mass)
        accarray = np.append(accarray,acceleration)
        
    # Create plots
    fig, ax1 = plt.subplots()
    
    ax2 = ax1.twinx()
    ax1.plot(timearray, altarray,'g-')
    ax2.plot(timearray, velarray,'b-')

    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Altitude (m)')
    ax2.set_ylabel('Velocity (m/s)')    
    plt.show()  
    
    return altarray,timearray,velarray,accarray
    
    
if __name__ == "__main__":
    
    alt_init = 70000. # m
    mass = 1000. # kg
    V_init = 200. # m/s
    C_D_chute = 0.8 # -
    A_chute = 19.6 # m^2
    g_venus = 8.87 # m/s^2
    dt = 0.1 # s
    mass_shield = 100. # kg
    t_shield = 30.0 # s
    lat = 0 # deg
    sol = 0 # deg
    SF = 1.2
    alt_end = 55000.
    altarray,timearray,velarray,accarray = ParachuteTrajectory(alt_init,mass,V_init,C_D_chute,A_chute,g_venus,dt,mass_shield,t_shield,lat,sol,SF,alt_end)

    
    print "Time to desired end altitude",timearray[-1], "seconds."
    print "Maximum acceleration is",np.max(abs(accarray)),"m/s^2 or",np.max(abs(accarray))/9.81,"g's."
    
    