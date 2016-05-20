# -*- coding: utf-8 -*-
"""
Created on Wed May 18 21:10:50 2016

@author: Yuyang
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import pylab
from mayavi import mlab 

import wing as w

"""open airfoil data"""
f = open('15012.txt', 'r')

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
x0up = d[:np.size(columnslist)/4, 0]
y0up = d[:np.size(columnslist)/4, 1]

#since xup is from 1 to 0, need to flip it
xup = np.zeros(np.size(columnslist)/4)
yup = np.zeros(np.size(columnslist)/4)
for i in range(np.size(columnslist)/4):
    xup[np.size(columnslist)/4 - 1 - i] = x0up[i]
    yup[np.size(columnslist)/4 - 1 - i] = y0up[i]
#xup = xup.tolist()
#yup = yup.tolist()
    
xlow = d[np.size(columnslist)/4:, 0]
ylow = d[np.size(columnslist)/4:, 1]
#xlow = xlow.tolist()
#ylow = ylow.tolist()

#data combined lists
#x0 = d[:, 0]
#y0 = d[:, 1]
x0 = np.concatenate((xup, xlow))
y0 = np.concatenate((yup, ylow))
#x0 = x0.tolist()
#y0 = y0.tolist()

"""Definitions"""
#dist boom to to y = 0, boom areas, neutral axis
def Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar):
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
    
    #init area of left and right cell
    Al = Ar = 0.
    
    #from 0.15x to 0.75x. start 9th point to 21th point. so probably 8 booms up
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
        #already positive
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
        #at spar, which is nr 3. boom
        if i == pos_spar - 1:
            Bup[i] = Bsup[i] + 2*Aflange
            Blow[i] = Bslow[i] + 2*Aflange
        else:
            Bup[i] = Bsup[i] + Aflange
            Blow[i] = Bslow[i] + Aflange
        
        #to calc neutral axis y coor
        ynaA.append(ytop[i]*Bup[i])
        ynaA.append(-ybot[i]*Blow[i])
        
        #calc left and right cell area
        if i < pos_spar - 1:
            Al = Al + (ytop[i] + ytop[i + 1] + ybot[i] + ybot[i + 1])/2.*pitch*sc[pos]**2
        elif i < nboom/2 - 1:
            Ar = Ar + (ytop[i] + ytop[i + 1] + ybot[i] + ybot[i + 1])/2.*pitch*sc[pos]**2
        
    #neutral axis, x axis is in the middle since pitch same for every boom
    xna = (x0[stop] - x0[start])/2 + x0[start]
    yna = sum(ynaA)/(sum(Bup) + sum(Blow))
    
    return ytop*sc[pos], ybot*sc[pos], Bup*sc[pos]**2, Blow*sc[pos]**2, \
            xna*sc[pos], yna, bup*sc[pos], blow*sc[pos], Al, Ar
    
#MOI
def MOI(nboom, tskin, pitch_coor, Aflange):
    #init lists
    Ixxlist = []
    Iyylist = []
    Ixylist = []
    
    for i in range(nboom/2):
        #get lists from Geometry
        ytop, ybot, Bup, Blow, xna, yna, bup, blow = Geometry(nboom, tskin, pitch_coor, \
                                            Aflange, pos, start, stop, pos_spar)[:-2]
        
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
def NormalStress(pos, i, start, pitch, sc):
    #get values from MOI
    Ixx, Iyy, Ixy = MOI(nboom, tskin, pitch_coor, Aflange)[:]
    
    #get value from Geometry
    xna, yna = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[4:6]
    
    Mx = (w.loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[2])[pos]
    My = (w.loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[3])[pos]
    
    #point of interest
    #x = (xup[start + i*pitch_coor] - xna)*sc[pos]
    x = (xup[i] - xna)*sc[pos]
    x_dlower = (xlow[i] - xna)*sc[pos]
    #y = (yup[start + i*pitch_coor] - yna)*sc[pos]
    y = (yup[i] - yna)*sc[pos]
    y_dlower = (ylow[i] - yna)*sc[pos]
    
    sigz = (Ixx*My - Ixy*Mx)/(Ixx*Iyy - Ixy**2)*x + (Iyy*Mx - Ixy*My)/ \
            (Ixx*Iyy - Ixy**2)*y
    sigzbot = (Ixx*My - Ixy*Mx)/(Ixx*Iyy - Ixy**2)*x_dlower + \
            (Iyy*Mx - Ixy*My)/(Ixx*Iyy - Ixy**2)*y_dlower
    return sigz, x, y, sigzbot, x_dlower, y_dlower

#open section shear
def openShearflow(sx, sy, Ixx, Iyy, Ixy, Bup, Blow, pitch_coor):
    #assume the skin takes shear only, so  tskin = 0
    qbuplist = []
    qblowlist = []
    for j in range(nboom/2):
        Bxup = Bup[j]*xup[start + j*pitch_coor]
        Byup = Bup[j]*yup[start + j*pitch_coor]
        qbup = -(Ixx*sx - Ixy*sy)/(Ixx*Iyy - Ixy**2)*Bxup - \
                (Iyy*sy - Ixy*sx)/(Ixx*Iyy - Ixy**2)*Byup
        qbuplist.append(qbup)
        
        Bxlow = Blow[j]*xlow[start + j*pitch_coor]
        Bylow = Blow[j]*ylow[start + j*pitch_coor]
        qblow = -(Ixx*sx - Ixy*sy)/(Ixx*Iyy - Ixy**2)*Bxlow - \
                (Iyy*sy - Ixy*sx)/(Ixx*Iyy - Ixy**2)*Bylow
        qblowlist.append(qblow)
    
    #cut 1-3 and 3-8, numbering from top left clockwise
    #assume shear flows positive counterclkwise
    #qbtot = sum(qbuplist) + sum(qblowlist)
    
    #shear flow at each section
    q0301 = 0.
    q0116 = qbuplist[0]
    q1614 = qblowlist[0]
    q1403 = qblowlist[2]
    q1409 = qblowlist[2]
    q0908 = qblowlist[7]
    q0803 = 0.
    
    return q0301, q0116, q1614, q1403, q1409, q0908, q0803

#rate of twist to calc qs   
def rateoftwist(G, t, sx, sy):
    ytop, ybot, Bup, Blow = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[:4]
    #get distance between booms 
    bup, blow = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[6:8]
    #get area of left and right cell
    Al, Ar = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[8:10]
    
    Ixx, Iyy, Ixy = MOI(nboom, tskin, pitch_coor, Aflange)[:]
    
    #distance between booms
    #ybot is positive
    L0116 = ytop[0] + ybot[0]
    L1614 = sum(blow[0:2])
    L1403 = ytop[2] - ybot[2]
    L0301 = sum(bup[0:2])
    L0314 = L1403
    L1409 = sum(blow[2:7])
    L0908 = ytop[7] + ybot[7]
    L0803 = sum(bup[2:7])
    
    q0301, q0116, q1614, q1403, q1409, q0908, q0803 = \
    openShearflow(sx, sy, Ixx, Iyy, Ixy, Bup, Blow, pitch_coor)[:] 
    
    #contribution by qb of left cell
    ctbtl = 1./(2*Al*G)*(q0116*L0116/t + q1614*L1614/t + q1403*L1403/t + q0301*L0301/t)
    #contribution by qb of right cell    
    ctbtr = 1./(2*Ar*G)*(q1403*L0314/t + q1409*L1409/t + q0908*L0908/t + q0803*L0803/t)
    
    #fill ctbtl and ctbtr back in dtheta/dz
#    dthetadzl = 1./(2*Al)*(qs1*(L0116/t + L1613/t + L1304/t + L0401/t) - qs2*L1304/t) + ctbtl
#    dthetadzr = 1./(2*Ar)*(-qs1*L1304/t + qs2*(L1304/t + L1309/t + L0908/t + L0804/t)) + ctbtr
    
    #to solve the quadratice eq of qs1 and qs2
#    (1./(2*Al)*(L0116/t + L1614/t + L1403/t + L0301/t) + 1./(2*Ar)*L1403/t)*qs1 + \
#    (-1./(2*Ar)*(L1403/t + L1409/t + L0908/t + L0803/t) - 1./(2*Al)*L1403/t)*qs2 \
#    = ctbtr - ctbtl
    
    #coefficient of the var for quad eq
    eq00 = (1./(2*Al)*(L0116/t + L1614/t + L1403/t + L0301/t) + 1./(2*Ar)*L1403/t)
    eq01 = (-1./(2*Ar)*(L1403/t + L1409/t + L0908/t + L0803/t) - 1./(2*Al)*L1403/t)
    eqr0 = ctbtr - ctbtl
    
    return L0116, L1614, L1403, L0301, L0314, L1409, L0908, L0803, ctbtl, ctbtr, eq00, eq01, eqr0
    
#torsion equilibruim to solve quadratic eq of qs1 and qs2
def torsionequivalence(sx, sy, pos, G, t):
    ybot = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[1]   
    
    #get boom area
    Bup, Blow = \
    Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[2:4]  
    
    yna = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[5]
    
    Al, Ar = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[8:10]
    
    #get boom moi info
    Ixx, Iyy, Ixy = MOI(nboom, tskin, pitch_coor, Aflange)[:]
    
    #get all boom distances
    L0116, L1614, L1403, L0301, L0314, L1409, L0908, L0803 = rateoftwist(G, t, sx, sy)[:8]

    #get qb of everywhere
    q0301, q0116, q1614, q1403, q1409, q0908, q0803 = \
    openShearflow(sx, sy, Ixx, Iyy, Ixy, Bup, Blow, pitch_coor)[:]   
    
    #moment equilibruim at middle spar
#    q0116*L0116*(xup[12] - xup[8])*sc[pos] + q1614*L1614*(ybot[1] - yna) + \
#    2*Al*qs1 + 2*Ar*qs2 = 0
    
    #2nd quadratic eq, derived from above
    #2*Al*qs1 + 2*Ar*qs2 = -q0116*L0116*(xup[12] - xup[8])*sc[pos] - q1614*L1614*(ybot[1] - yna)
    
    #coefficient from var in 2nd quadratic eq
    eq10 = 2*Al
    eq11 = 2*Ar
    eqr1 = -q0116*L0116*(xup[12] - xup[8])*sc[pos] - q1614*L1614*(ybot[1] - yna)
    
    return eq10, eq11, eqr1  
        
"""Parameters. Change this if needed!"""
#total nr of booms
nboom = 16
#start at the 9th element
start = 8
#stop at the 23th element
stop = 22
#position/section along span
pos = 1
#the nth boom there's a spar, 3 cuz it's at 0.25c
pos_spar = 3

tskin = 0.005
t = tskin
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

sx = 0.
sy = 1000
G = 41.1*10**9

#scaling of airfoil
sc = w.loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[-1]

"""Excution"""
#to get valus of sigz
#init sigz list
#sigzlist = np.zeros(np.size(sc) - 1)
sigzlist = [[]]
#make the list 3 columns
sigzlist.append([])
sigzlist.append([])

#list for x and y and pos and sigz, upper part
xylist = [[]]
xylist.append([])
xylist.append([])
xylist.append([])

#list for x and y and pos and sigz, lower part
xybotlist = [[]]
xybotlist.append([])
xybotlist.append([])
xybotlist.append([])

#init z for the color map
z = pylab.zeros([np.size(sc) - 1, nboom/2])

#init list of q
qlist = [[]]
for i in range(6):
    qlist.append([])

for pos in range(np.size(sc) - 1):
    #distance from boom to y = 0
    ytop = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[0]
    #this is negative already
    ybot = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[1]
    
    Bup = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[2]
    Blow = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[3]
    
    #moment of inertia
    Ixx = MOI(nboom, tskin, pitch_coor, Aflange)[0]
    Iyy = MOI(nboom, tskin, pitch_coor, Aflange)[1]
    Ixy = MOI(nboom, tskin, pitch_coor, Aflange)[2]
    
    #normal stress
    for i in range(8):        
        z[pos, i] = NormalStress(pos, start + i*pitch_coor, start, pitch_coor, sc)[0]
    
#    #shear flow calc of each section
#    eq00, eq01, eqr0 = rateoftwist(G, t, sx, sy)[-3:]
#    eq10, eq11, eqr1 = torsionequivalence(sx, sy, pos, G, t)[:]
#    
#    m = np.array([[eq00, eq01],[eq10, eq11]])
#    n = np.array([eqr0, eqr1])
#    qs1, qs2 = np.linalg.solve(m, n)[:]
#    print qs1, qs2
#    
#    q0301, q0116, q1614 = \
#    openShearflow(sx, sy, Ixx, Iyy, Ixy, Bup, Blow, pitch_coor)[:3] - qs1
#    
#    qlist[0].append(q0301)
#    qlist[1].append(q0116)
#    qlist[2].append(q1614)    
#    
#    q1403 = openShearflow(sx, sy, Ixx, Iyy, Ixy, Bup, Blow, pitch_coor)[3] + qs1 - qs2
#    
#    qlist[3].append(q1403)
#    
#    q1409, q0908, q0803 = \
#    openShearflow(sx, sy, Ixx, Iyy, Ixy, Bup, Blow, pitch_coor)[4:] + qs2
#    
#    qlist[4].append(q1409)
#    qlist[5].append(q0908)
#    qlist[6].append(q0803)        
    
    for i in range(np.size(xup)):
        #upper part
        x = NormalStress(pos, i, start, pitch_coor, sc)[1]
        y = NormalStress(pos, i, start, pitch_coor, sc)[2]
        xylist[0].append(x*10.)
        xylist[1].append(y*10.)
        xylist[2].append(pos)
        sigz = NormalStress(pos, i, start, pitch_coor, sc)[0]
        xylist[3].append(-sigz)
        
        #add lower part together
        x_dlower = NormalStress(pos, i, start, pitch_coor, sc)[4]
        y_dlower = NormalStress(pos, i, start, pitch_coor, sc)[-1]
        xylist[0].append(x_dlower*10.)
        xylist[1].append(y_dlower*10.)
        xylist[2].append(pos)
        sigzbot = NormalStress(pos, i, start, pitch_coor, sc)[3]
        xylist[3].append(-sigzbot)
        
"""Plot 2d wingbox"""
#color map plot
pos = np.linspace(0, np.size(sc) - 2, np.size(sc) - 1)
i = np.linspace(0, 7, 8)
#gotta flip i and pos for some reason
sigmap = pylab.pcolor(i, pos, z)

pylab.colorbar(sigmap)
sigmap.colorbar.set_label('Normal stress [Pa]')
pylab.xlabel('Boom sections')
pylab.ylabel('Span sections root to tip')
pylab.show()

##3d plot of wing
#def plot_mayavi(data):
#    fig = mlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0))
#    obj = mlab.points3d(data[0], data[1], data[2], data[3], line_width = 1.0, mode="cube", scale_mode="none",scale_factor=0.5)
#    bar = mlab.colorbar(object=obj, title="Stress (in Pa)", orientation="vertical")
#    mlab.view(50,50)
#    mlab.roll(5)
#    #mlab.savefig(os.path.dirname(os.path.realpath(__file__))+'\\sim_outputs\\plot_mayavi.jpg', size=[5000,3000])
#    mlab.show()  
#    
#plot_mayavi(xylist)   
#
#"""plot wingbox geometry"""
##preset pos for plotting purpose
#pos = 1
#
##distance from boom to y = 0
#ytop = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[0]
##this is negative already
#ybot = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[1]
#
##neutral axis
#xna = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[4]
#yna = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[5]
#
##init x1 list, which is x axis
#x1 = np.zeros(nboom/2)
#for i in range(nboom/2):
#    x1[i] = x0[start + i*pitch_coor]*sc[pos]
#
#plt.plot(x1, ytop, linewidth = 3, linestyle = '-', c = 'r', label = 'Wingbox')
#plt.plot(x1, -ybot, linewidth = 3, linestyle = '-', c = 'r')
#
##plot neutral axis
#plt.plot(np.ones(nboom/2)*xna, np.linspace(-0.3, 0.3, nboom/2), 'b--')
#plt.plot(np.linspace(-0.2, 1.4, nboom/2), np.ones(nboom/2)*yna, 'b--')
#
##plot start and end spars
#plt.plot(np.ones(nboom/2)*x1[0], np.linspace(-ybot[0], ytop[0], nboom/2), \
#        linewidth = 3, linestyle = '-', c = 'r')
#plt.plot(np.ones(nboom/2)*x1[-1], np.linspace(-ybot[-1], ytop[-1], nboom/2), \
#        linewidth = 3, linestyle = '-', c = 'r')
#        
##plot middle spar which is at 4th boom from left
#plt.plot(np.ones(nboom/2)*x1[2], np.linspace(-ybot[2], ytop[2], nboom/2), \
#        linewidth = 3, linestyle = '-', c = 'r')
#
##plot airfoil
#plt.plot(sc[pos]*xup, sc[pos]*yup, 'g-')
#plt.plot(sc[pos]*xlow, sc[pos]*ylow, 'g-')
#plt.axis('equal')
#
##plot booms
#i = 0
#while i < nboom/2:     
#    circletop=plt.Circle((x1[i], ytop[i]), .01,color='r', clip_on = False)
#    fig = plt.gcf()
#    fig.gca().add_artist(circletop)
#    circletop=plt.Circle((x1[i], -ybot[i]), .01,color='r', clip_on = False)
#    fig = plt.gcf()
#    fig.gca().add_artist(circletop)
#    i = i + 1
#
##plt.ylim(-0.5, 0.5)
#plt.xlabel('Chord length [m]')
#plt.ylabel('Camber [m]')
#plt.legend(loc = 1)
#plt.tight_layout()
#plt.show()

##plot stress curve
#plt.plot(np.linspace(0, pos, pos + 1), sigzlist/10**6)
#plt.xlabel('Position on the wing span [-]')
#plt.ylabel('Normal stress [MPa]')
#plt.show()

##plot shear flow diagram
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
##0116
#X1 = np.linspace(xup[start], xlow[start], 10)
#Y1 = np.linspace(yup[start], ylow[start], 10)
#Z1 = np.ones(10)*qlist[1][0]
#ax.plot_wireframe(X1, Y1, Z1, color = 'g')
#
##1614
#X2 = np.linspace(xlow[start], xlow[start + pitch_coor*(pos_spar - 1)], 10)
#Y2 = np.linspace(ylow[start], ylow[start + pitch_coor*(pos_spar - 1)], 10)
#"""This is wrong, the sign is wrong"""
#Z2 = -np.ones(10)*qlist[2][0]
#ax.plot_wireframe(X2, Y2, Z2, color = 'g')
#
##1403
#X3 = np.linspace(xup[start + pitch_coor*(pos_spar - 1)], xlow[start + pitch_coor*(pos_spar - 1)], 10)
#Y3 = np.linspace(yup[start + pitch_coor*(pos_spar - 1)], ylow[start + pitch_coor*(pos_spar - 1)], 10)
#Z3 = np.ones(10)*qlist[3][0]
#ax.plot_wireframe(X3, Y3, Z3, color = 'r')
#
##0301
#X4 = np.linspace(xup[start], xup[start + pitch_coor*(pos_spar - 1)], 10)
#Y4 = np.linspace(yup[start], yup[start + pitch_coor*(pos_spar - 1)], 10)
#Z4 = np.ones(10)*qlist[0][0]
#ax.plot_wireframe(X4, Y4, Z4, color = 'g')
#
##1409
#X5 = np.linspace(xlow[start + pitch_coor*(pos_spar - 1)], xlow[stop], 10)
#Y5 = np.linspace(ylow[start + pitch_coor*(pos_spar - 1)], ylow[stop], 10)
#Z5 = np.ones(10)*qlist[4][0]
#ax.plot_wireframe(X5, Y5, Z5, color = 'r')
#
##0908
#X6 = np.linspace(xup[stop], xlow[stop], 10)
#Y6 = np.linspace(yup[stop], ylow[stop], 10)
#Z6 = np.ones(10)*qlist[5][0]
#ax.plot_wireframe(X6, Y6, Z6, color = 'r')
#
##0803
#X7 = np.linspace(xup[start + pitch_coor*(pos_spar - 1)], xup[stop], 10)
#Y7 = np.linspace(yup[start + pitch_coor*(pos_spar - 1)], yup[stop], 10)
#Z7 = np.ones(10)*qlist[6][0]
#ax.plot_wireframe(X7, Y7, Z7, color = 'r')
#
#"""base geometry projection"""
#X8 = np.linspace(xup[start], xlow[start], 10)
#Y8 = np.linspace(yup[start], ylow[start], 10)
#Z8 = np.zeros(10)
#ax.plot_wireframe(X8, Y8, Z8, color = 'b', linestyles = '--')
#
##1614
#X9 = np.linspace(xlow[start], xlow[start + pitch_coor*(pos_spar - 1)], 10)
#Y9 = np.linspace(ylow[start], ylow[start + pitch_coor*(pos_spar - 1)], 10)
#Z9 = np.zeros(10)
#ax.plot_wireframe(X9, Y9, Z9, color = 'b', linestyles = '--')
#
##1403
#X10 = np.linspace(xup[start + pitch_coor*(pos_spar - 1)], xlow[start + pitch_coor*(pos_spar - 1)], 10)
#Y10 = np.linspace(yup[start + pitch_coor*(pos_spar - 1)], ylow[start + pitch_coor*(pos_spar - 1)], 10)
#Z10 = np.zeros(10)
#ax.plot_wireframe(X10, Y10, Z10, color = 'b', linestyles = '--')
#
##0301
#X11 = np.linspace(xup[start], xup[start + pitch_coor*(pos_spar - 1)], 10)
#Y11 = np.linspace(yup[start], yup[start + pitch_coor*(pos_spar - 1)], 10)
#Z11 = np.zeros(10)
#ax.plot_wireframe(X11, Y11, Z11, color = 'b', linestyles = '--')
#
##1409
#X12 = np.linspace(xlow[start + pitch_coor*(pos_spar - 1)], xlow[stop], 10)
#Y12 = np.linspace(ylow[start + pitch_coor*(pos_spar - 1)], ylow[stop], 10)
#Z12 = np.zeros(10)
#ax.plot_wireframe(X12, Y12, Z12, color = 'b', linestyles = '--')
#
##0908
#X13 = np.linspace(xup[stop], xlow[stop], 10)
#Y13 = np.linspace(yup[stop], ylow[stop], 10)
#Z13 = np.zeros(10)
#ax.plot_wireframe(X13, Y13, Z13, color = 'b', linestyles = '--')
#
##0803
#X14 = np.linspace(xup[start + pitch_coor*(pos_spar - 1)], xup[stop], 10)
#Y14 = np.linspace(yup[start + pitch_coor*(pos_spar - 1)], yup[stop], 10)
#Z14 = np.zeros(10)
#ax.plot_wireframe(X14, Y14, Z14, color = 'b', linestyles = '--')
#
#ax.set_xlabel('Chord location [-]')
#ax.set_ylabel('Camber location [-]')
#ax.set_zlabel('Shear flow [N/m]')

plt.tight_layout()

        
        
        
    
    
    
