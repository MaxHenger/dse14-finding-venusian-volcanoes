# -*- coding: utf-8 -*-
"""
Created on Thu May 19 17:09:00 2016

@author: MaxHenger
"""

import Atmosphere as atmosphere
#import Gravity as gravity

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.integrate as scpint

def ProbeAtmosphere():
    atm = atmosphere.Atmosphere()
    
    z = np.linspace(0.0, 70000.0, 200)
    rho = atm.density(z, 0, 0)
    vel = atm.velocityZonal(z, 0, 0)
    
    deltaV = [-30, -20, -10, 0, 10, 20, 30]
    
    fig = plt.figure()
    
    curMax = -1e9
    curMin = 1e9
    axes = []
    for i in range(0, len(deltaV)):
        curQ = 0.5 * rho[1] * np.power(vel[1] + deltaV[i], 2.0)
        ax = fig.add_subplot(100 + 10 * len(deltaV) + i + 1)
        ax.plot(curQ, z / 1000)
        ax.grid(True)
        ax.set_title('deltaV = ' + str(deltaV[i]))
        axes.append(ax)
        if max(curQ) > curMax:
            curMax = max(curQ)
            
        if min(curQ) < curMin:
            curMin = min(curQ)
        
    for i in range(0, len(axes)):
        axes[i].set_xlim(curMin, curMax)
        
    fig.tight_layout()
    print('v at 30 =', atm.velocityZonal(30000, 0, 0))
    print('v at 40 =', atm.velocityZonal(40000, 0, 0))
    print('v at 60 =', atm.velocityZonal(60000, 0, 0))
    print('rho at 30 =', atm.density(30000, 0, 0))
    print('rho at 60 =', atm.density(40000, 0, 0))
    print('rho at 60 =', atm.density(60000, 0, 0))
    
#ProbeAtmosphere()

def ProbeGravity():
    z = np.linspace(0.0, 70000, 500)
    #grav = gravity.Gravity()
    #cur = grav(z, 0, 0)
    
    plt.plot(z, cur)
    
#ProbeGravity()
    
def PlotAtmosphere():
    atm = atmosphere.Atmosphere()
    
    z = np.linspace(0.0, 70000, 500)
    rho = atm.density(z, 0, 0)
    vel = atm.velocityZonal(z, 0, 0)
    
    fig = plt.figure()
    
    # Draw pressure and density with error bounds
    axVel = fig.add_axes([0.08, 0.10, 0.30, 0.80])
    axDen = axVel.twiny()
    
    axVel.plot(vel[0], z / 1000, 'g--')
    lineVel, = axVel.plot(vel[1], z / 1000, 'g')
    axVel.plot(vel[2], z / 1000, 'g--')
    
    axVel.set_xlabel(r'$V\;[m/s]$', fontsize=15, color='g')
    axVel.set_ylabel(r'$h\;[km]$', fontsize=15)
    axVel.tick_params(axis='x', colors='g')
    axVel.grid(True)
    
    axDen.plot(np.log10(rho[0]), z / 1000, 'b--')
    lineDen, = axDen.plot(np.log10(rho[1]), z / 1000, 'b')
    axDen.plot(np.log10(rho[2]), z / 1000, 'b--')
    axDen.set_xlabel(r'$\mathrm{log}(p)\;[Pa]$', fontsize=15, color='b')
    axDen.tick_params(axis='x', colors='b')
    
    # Draw dynamic pressure as a function of delta V
    deltaV = np.linspace(-45, 45, 250)

    meanDynPres = np.zeros([len(z), len(deltaV)])
    
    for iZ in range(0, len(z)):
        for iDV in range(0, len(deltaV)):
            meanDynPres[iZ, iDV] = np.log10(0.5 * rho[1][iZ] * (vel[1][iZ] + deltaV[iDV])**2.0)
            
    axDyn = fig.add_axes([0.48, 0.10, 0.38, 0.80])
    gnuplot2 = plt.get_cmap('gnuplot2')
    
    vmin = np.min(meanDynPres)
    vmax = np.max(meanDynPres)
    norm = mpl.colors.Normalize(0.0, vmax)
    axDyn.imshow(meanDynPres, extent=[deltaV[0], deltaV[-1], z[0] / 1000, z[-1] / 1000],
                origin='lower', aspect='auto', cmap=gnuplot2, norm=norm)
    contourLower = axDyn.contour(deltaV, z / 1000, meanDynPres, [1.5], colors='w')
    axDyn.clabel(contourLower, fontsize=9)
    contourUpper = axDyn.contour(deltaV, z / 1000, meanDynPres, np.linspace(2.5, 5.0, 6), colors='k')
    axDyn.clabel(contourUpper, fontsize=9)
    axDyn.set_xlabel(r'$\Delta V \; [m/s]$', fontsize=15)
    axDyn.set_ylabel(r'$h\;[km]$', fontsize=15)
    axDyn.grid(True)
    
    # Plot the colormap
    axMap = fig.add_axes([0.88, 0.10, 0.03, 0.80])
    cbb = mpl.colorbar.ColorbarBase(axMap, cmap=gnuplot2, norm=norm)
    cbb.set_label(r'$\mathrm{log}(q_\infty)\;[Pa]$', rotation=90, fontsize=15)
    
PlotAtmosphere()

def PlotWingArea():
    #Cr = np.linspace(1.0, 5.0, 100)
    #CtOverCr = np.linspace(0.1, 1.0, 100)
    q = 1000
    Cr = np.linspace(0.5, 4.5, 10)
    CtOverCr = np.linspace(0.1, 1.0, 10)
    Area = np.zeros([len(CtOverCr), len(Cr)])
    Cl = np.zeros([len(CtOverCr), len(Cr)])
    ClTerm = np.zeros([len(CtOverCr), len(Cr)])
    W = 750 * 8.7
    R = 2.5
    bottom = 0.13
    top = 0.93
    
    for iCr in range(0, len(Cr)):
        for iCtOverCr in range(0, len(CtOverCr)):
            Ct = Cr[iCr] * CtOverCr[iCtOverCr]
            
            # Solve for wingspan
            b1Overb2 = np.linspace(0.5, 1.5, 100)
            recalculatedb1Overb2 = np.sqrt(R**2.0 - 0.25 * Ct**2.0) / \
                np.sqrt(R**2.0 - 0.25 * np.power(Cr[iCr] + (Ct - Cr[iCr]) / (b1Overb2 + 2.0), 2.0))
                
            iBest = np.argmin(abs(b1Overb2 - recalculatedb1Overb2))
            b1 = np.sqrt(R**2.0 - 0.25 * Ct**2.0)
            b2 = np.sqrt(R**2.0 - 0.25 * (Cr[iCr] + (Ct - Cr[iCr]) / (b1Overb2[iBest] + 2.0))**2.0)
            b = 2 * (b1 + 2 * b2)
            
            # calculate area and lift parameter
            x = np.linspace(0, b / 2, 200)
            yArea = Cr[iCr] + (Ct - Cr[iCr]) / (b / 2) * x
            yLift = yArea * np.sqrt(1 - np.power(2 * x / b, 2.0))
            
            Area[iCtOverCr, iCr] = 2.0 * scpint.simps(yArea, x)
            ClTerm[iCtOverCr, iCr] = scpint.simps(yLift, x)
            Cl[iCtOverCr, iCr] = W / (2.0 * q * ClTerm[iCtOverCr, iCr])
            
            
    gnuplot2 = plt.get_cmap("gnuplot2")
    
    # Plot the wing area
    fig = plt.figure()
    axArea = fig.add_axes([0.10, bottom, 0.32, top - bottom])
    
    vmin = np.min(Area)
    vmax = np.max(Area)
    norm = mpl.colors.Normalize(vmin, vmax)
    axArea.imshow(Area, extent=[Cr[0], Cr[-1], CtOverCr[0], CtOverCr[-1]],
                  origin='lower', aspect='auto', cmap=gnuplot2, norm=norm)
    axArea.set_xlabel(r'$c_r\;[m]$', fontsize=14)
    axArea.set_ylabel(r'$\lambda\;=\;c_t/c_r$', fontsize=14)
    contourAreaLower = axArea.contour(Cr, CtOverCr, Area, np.linspace(10, 20, 2), colors='w')
    axArea.clabel(contourAreaLower, fontsize=9)
    contourAreaUpper = axArea.contour(Cr, CtOverCr, Area, np.linspace(25, 60, 8), colors='k')
    axArea.clabel(contourAreaUpper, fontsize=9)
    
    axAreaMap = fig.add_axes([0.44, bottom, 0.03, top - bottom])
    cbAreaMap = mpl.colorbar.ColorbarBase(axAreaMap, cmap=gnuplot2, norm=norm)
    cbAreaMap.set_label(r'$S\;[m^2]$', rotation=90, fontsize=15)
    
    # Plot the design lift coefficient
    axCl = fig.add_axes([0.60, bottom, 0.28, top - bottom])
    
    vmin = np.min(Cl)
    vmax = np.max(Cl)
    norm = mpl.colors.Normalize(vmin, vmax)
    axCl.imshow(Cl, extent=[Cr[0], Cr[-1], CtOverCr[0], CtOverCr[-1]],
                origin='lower', aspect='auto', cmap=gnuplot2, norm=norm)
    axCl.set_xlabel(r'$c_r\;[m]$', fontsize=14)
    axCl.set_ylabel(r'$\lambda\;=\;c_t/c_r$', fontsize=14)
    
    contoursLower = np.linspace(0, 0.04, 9)
    contoursUpper = np.linspace(0.05, 0.20, 4)
    
    if q < 5000:
        contoursLower = np.linspace(0, 0.4, 9)
        contoursUpper = np.linspace(0.5, 1.0, 3)
        
    contourClLower = axCl.contour(Cr, CtOverCr, Cl, contoursLower, colors='w')
    axCl.clabel(contourClLower, fontsize=9)
    contourClUpper = axCl.contour(Cr, CtOverCr, Cl, contoursUpper, colors='w')
    axCl.clabel(contourClUpper, fontsize=10)
    
    axClMap = fig.add_axes([0.89, bottom, 0.03, top - bottom])
    cbClMap = mpl.colorbar.ColorbarBase(axClMap, cmap=gnuplot2, norm=norm)
    cbClMap.set_label(r'$c_l$', rotation=90, fontsize=15)
    
PlotWingArea()

def GetWingspan(Cr, Ct, R):
    # Solve for wingspan
    b1Overb2 = np.linspace(0.5, 1.5, 1000)
    recalculatedb1Overb2 = np.sqrt(R**2.0 - 0.25 * Ct**2.0) / \
        np.sqrt(R**2.0 - 0.25 * np.power(Cr + (Ct - Cr) / (b1Overb2 + 2.0), 2.0))
        
    iBest = np.argmin(abs(b1Overb2 - recalculatedb1Overb2))
    b1 = np.sqrt(R**2.0 - 0.25 * Ct**2.0)
    b2 = np.sqrt(R**2.0 - 0.25 * (Cr + (Ct - Cr) / (b1Overb2[iBest] + 2.0))**2.0)
    b = b1 + 2 * b2
    
    print('b1 =', b1)
    print('b2 =', b2)
    print('b  =', b)
    
#GetWingspan(3.5, 1.5, 2.5)
    
def ProbeWingSizes():
    Cr = np.linspace(0.5, 3.5, 100)
    CtCr = np.linspace(0.1, 1.0, 100)
    b = 15
    q_inf = [500, 6000, 8000]
    Cl = np.zeros([len(Cr), len(CtCr), len(q_inf)])
    W = 750*8.8
    
    for iCr in range(0, len(Cr)):
        for iRatio in range(0, len(CtCr)):
            for iQ in range(0, len(q_inf)):
                Ct = Cr[iCr] * CtCr[iRatio]
                
                x = np.linspace(0, b / 2, 200)
                y = Cr[iCr] + (Ct - Cr[iCr]) / (b / 2) * x * np.sqrt(1 - np.power(2 * x / b, 2))
                integral = scpint.simps(y, x)
                Cl[iCr, iRatio, iQ] = W / (2.0 * q_inf[iQ] * integral)
    
    fig = plt.figure()
    
    cmap = plt.get_cmap('jet')
    
    for i in range(0, len(q_inf)):
        ax = fig.add_subplot(100 + 10 * len(q_inf) + i + 1)
        ax.imshow(Cl[:, :, i], extent=[Cr[0], Cr[-1], CtCr[0], CtCr[-1]], origin='lower', aspect='auto', cmap=cmap)
        ax.set_title('q_inf = ' + str(q_inf[i]))
        
    #fig.colorbar(cmap)
        
#ProbeWingSizes()