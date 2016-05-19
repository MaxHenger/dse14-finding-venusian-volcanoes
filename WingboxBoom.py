# -*- coding: utf-8 -*-
"""
Created on Wed May 18 21:10:50 2016

@author: Yuyang
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches
import wing as w
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

#open airfoil data
f = open('airfoil.txt', 'r')

c = [0,0]

columnslist = []

for lines in f:
    lines = lines.strip()
    columns = lines.split()
    columnslist.append(columns)
for i in range(np.size(columnslist)/2):
    a = [float(j) for j in columnslist[i]]
    c = np.vstack((c,a))

d = c[1:,:]
print ''
print 'Airfoil dataset: '
print d

#data lists by top and bottom
xup = d[:np.size(columnslist)/4, 0]
yup = d[:np.size(columnslist)/4, 1]
#xup = xup.tolist()
#yup = yup.tolist()
xlow = d[np.size(columnslist)/4:, 0]
ylow = d[np.size(columnslist)/4:, 1]
#xlow = xlow.tolist()
#ylow = ylow.tolist()

#data combined lists
x0 = d[:, 0]
y0 = d[:, 1]
#x0 = x0.tolist()
#y0 = y0.tolist()

#dist boom to to y = 0, boom areas, neutral axis
def Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop):
    #assume top boom nr = bottom
    #init lists
    bup = np.zeros(nboom/2)
    blow = np.zeros(nboom/2)
    ytop = np.zeros(nboom/2)
    ybot = np.zeros(nboom/2)
    Bsup = np.zeros(nboom/2)
    Bslow = np.zeros(nboom/2)
    Bup = np.zeros(nboom/2)
    Blow = np.zeros(nboom/2)
    ynaA = []
    
    #from 0.15x to 0.75x. startth point to 21th point. so probably 8 booms up
    pitch = (x0[20] - x0[8])/(nboom/2 - 1)
    for i in range(nboom/2):
        #direct distance between booms, use triangle rule. 10 cuz start at 10th pt
        if i < nboom/2 - 1: 
            bup[i] = np.sqrt((yup[start + (i + 1)*pitch_coor] - yup[start + i*pitch_coor])**2 + pitch**2)
            blow[i] = np.sqrt((ylow[start + (i + 1)*pitch_coor] - ylow[start + i*pitch_coor])**2 + pitch**2)        
        else:
            bup[i] = bup[i - 1]
            blow[i] = blow[i - 1]
        
        # distance to y = 0, start cuz start at 10th pt
        ytop[i] = yup[start + i*pitch_coor]
        ybot[i] = abs(ylow[start + i*pitch_coor])
        
        #boom virtual area + skin
        if i == 0:
            Bsup[i] = tskin*bup[i]/6.*(2 + ytop[i + 1]/ytop[i]) + tskin*(ytop[i] + abs(ybot[i]))/6*(2 - 1)
            Bslow[i] = tskin*blow[i]/6.*(2 + ybot[i + 1]/ybot[i]) + tskin*(ytop[i] + abs(ybot[i]))/6*(2 - 1)
        elif i == nboom/2 - 1:
            Bsup[i] = tskin*bup[i]/6.*(2 + ytop[i]/ytop[i - 1]) + tskin*(ytop[i] + abs(ybot[i]))/6*(2 - 1)
            Bslow[i] = tskin*blow[i]/6.*(2 + ybot[i]/ybot[i - 1]) + tskin*(ytop[i] + abs(ybot[i]))/6*(2 - 1)
        else:
            Bsup[i] = tskin*bup[i]/6.*(2 + ytop[i + 1]/ytop[i])
            Bslow[i] = tskin*blow[i]/6.*(2 + ybot[i + 1]/ybot[i])
        
        #boom total area
        if i == nboom/4 - 1:
            Bup[i] = Bsup[i] + 2*Aflange
            Blow[i] = Bslow[i] + 2*Aflange
        else:
            Bup[i] = Bsup[i] + Aflange
            Blow[i] = Bslow[i] + Aflange
        
        #to calc neutral axis y coor
        ynaA.append(ytop[i]*Bup[i])
        ynaA.append(-ybot[i]*Blow[i])
    
    #neutral axis, x axis is in the middle since pitch same for every boom
    xna = (x0[stop] - x0[start])/2 + x0[start]
    yna = sum(ynaA)/(sum(Bup) + sum(Blow))
            
    return ytop*sc[pos], ybot*sc[pos], Bup*sc[pos]**2, Blow*sc[pos]**2, xna*sc[pos], yna
    
#MOI
def MOI(nboom, tskin, pitch_coor, Aflange):
    #init lists
    Ixxlist = []
    Iyylist = []
    Ixylist = []
    
    for i in range(nboom/2):
        #get lists from Geometry
        ytop, ybot, Bup, Blow, xna, yna = Geometry(nboom, tskin, pitch_coor, \
                                            Aflange, pos, start, stop)[:]
#        ybot = Geometry(nboom, tskin, pitch_coor, Aflange)[1]
#        Bup = Geometry(nboom, tskin, pitch_coor, Aflange)[2]
#        Blow = Geometry(nboom, tskin, pitch_coor, Aflange)[3]
#        xna = Geometry(nboom, tskin, pitch_coor, Aflange)[-2]
#        yna = Geometry(nboom, tskin, pitch_coor, Aflange)[-1]
        
        #Ixx of each boom
        Ixxlist.append((ytop[i] - yna)**2*Bup[i])
        Ixxlist.append((-ybot[i] - yna)**2*Blow[i])
        
        #Iyy of each boom
        Iyylist.append((xup[i] - xna)**2*Bup[i])
        Iyylist.append((xlow[i] - xna)**2*Blow[i])
        
        #Ixy of each boom, not 0 if boom not symmetric
        Ixylist.append((x0[start + i*pitch_coor] - xna)*(ytop[i] - yna)*Bup[i])
        Ixylist.append((x0[start + i*pitch_coor] - xna)*(-ybot[i] - yna)*Blow[i])
        
    #total Ixx
    Ixx = sum(Ixxlist)
    Iyy = sum(Iyylist)
    Ixy = sum(Ixylist)
    
    return Ixx, Iyy, Ixy

#Normal stress
def NormalStress(pos):
    #get values from MOI
    Ixx = MOI(nboom, tskin, pitch_coor, Aflange)[0]
    Iyy = MOI(nboom, tskin, pitch_coor, Aflange)[1]
    Ixy = MOI(nboom, tskin, pitch_coor, Aflange)[-1]
    
    #get value from Geometry
    xna, yna = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop)[-2:]
    
    Mx = (w.loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[2])[pos]
    My = 0.
    
    #point of interest
    x = (xup[1] - xna)*sc[pos]
    y = (yup[1] - yna)*sc[pos]
    
    sigz = (Ixx*My - Ixy*Mx)/(Ixx*Iyy - Ixy**2)*x + (Iyy*Mx - Ixy*My)/ \
            (Ixx*Iyy - Ixy**2)*y
            
    return sigz

"""change this if needed!"""
#total nr of booms
nboom = 16
#start at the 9th element
start = 8
#stop at the 21th element
stop = 20
#position/section along span
pos = 1

tskin = 0.005
pitch_coor = 2
Aflange = 300*10**(-6)

#from wing.py
b=50.
Cl=10.
Cd=1.
rho=3.
V=40.       #velocity in m/s
q=1./2*rho*V**2
A=10.
taper=0.5
F=0.006
tw=0.003
g=8.
S=10.
rhowing=3000.

#scaling of airfoil
sc = w.loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[-1]

##to get valus of sigz
#for pos in range(np.size(sc) - 1):
#    Ixx = MOI(nboom, tskin, pitch_coor, Aflange)[0]
#    Ixy = MOI(nboom, tskin, pitch_coor, Aflange)[1]
#
#    sigz = NormalStress(pos)
#
#    print -sigz

ytop = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop)[0]
ybot = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop)[1]

Ixx = MOI(nboom, tskin, pitch_coor, Aflange)[0]
Ixy = MOI(nboom, tskin, pitch_coor, Aflange)[1]

sigz = NormalStress(pos)

#plot wingbox
#init x1 list
x1 = np.zeros(nboom/2)
for i in range(nboom/2):
    x1[i] = x0[start + i*pitch_coor]*sc[pos]

plt.plot(x1, ytop, linewidth = 3, linestyle = '-', c = 'r')
plt.plot(x1, -ybot, linewidth = 3, linestyle = '-', c = 'r')

#plot start and end spars
plt.plot(np.ones(nboom/2)*x1[0], np.linspace(-ybot[0], ytop[0], nboom/2), \
        linewidth = 3, linestyle = '-', c = 'r')
plt.plot(np.ones(nboom/2)*x1[-1], np.linspace(-ybot[-1], ytop[-1], nboom/2), \
        linewidth = 3, linestyle = '-', c = 'r')
        
#plot middle spar which is at 4th boom from left
plt.plot(np.ones(nboom/2)*x1[3], np.linspace(-ybot[3], ytop[3], nboom/2), \
        linewidth = 3, linestyle = '-', c = 'r')

#plot airfoil
plt.plot(sc[pos]*xup, sc[pos]*yup, 'g-')
plt.plot(sc[pos]*xlow, sc[pos]*ylow, 'g-')
plt.axis('equal')

#plot booms
i = 0
while i < nboom/2:     
    circletop=plt.Circle((x1[i], ytop[i]), .01,color='r', clip_on = False)
    fig = plt.gcf()
    fig.gca().add_artist(circletop)
    circletop=plt.Circle((x1[i], -ybot[i]), .01,color='r', clip_on = False)
    fig = plt.gcf()
    fig.gca().add_artist(circletop)
    i = i + 1

#plt.ylim(-0.5, 0.5)
plt.show()
        

        
        
        
    
    
    