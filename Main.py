# -*- coding: utf-8 -*-
"""
Created on Wed Jun 01 02:10:07 2016

@author: Yuyang
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches
from matplotlib import cm
import matplotlib as mpl
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import pylab
from mayavi import mlab 

#import WinboxBoom for bending and shear, import wing for wing loading and moment
import WingboxBoom as WB
import wing as w

"""Parameters. Change this if needed!"""
        
tskin = 0.005
t = tskin
Aflange = 300*10**(-6)

#from wing.py
b=50.
Cl=0.175
Cd=0.0075
rho=3.
V=40.       #velocity in m/s
q=1000.
A=10.
taper=0.55
F=0.
tw=0.005
g=8.
S=35.
rhowing=3000.

"""need to change this"""
sx = w.loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[1]
sy = w.loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[0]
G = 41.1*10**9

#scaling of airfoil
sc = w.loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[-1]

ver = '8/2'

#print WB.openShearflow(sx, sy, ver)

"""Excution"""
##to get valus of sigz
##init sigz list
##sigzlist = np.zeros(np.size(sc) - 1)
#sigzlist = [[]]
##make the list 3 columns
#sigzlist.append([])
#sigzlist.append([])
#
##list for x and y and pos and sigz, upper part
#xylist = [[]]
#xylist.append([])
#xylist.append([])
#xylist.append([])
#
##list for x and y and pos and sigz, lower part
#xybotlist = [[]]
#xybotlist.append([])
#xybotlist.append([])
#xybotlist.append([])
#
##init z(sigma) for the color map
#z = pylab.zeros([np.size(sc) - 1, nboom/2])
#zbot = pylab.zeros([np.size(sc) - 1, nboom/2])
#
##init list of q
#qlist = [[]]
#for i in range(nboom - 2 + nspar):
#    qlist.append([])
#

#get config from config def
start, stop, nboom, pitch, pitch_coor, pos_spar1, pos_spar2, nspar = WB.config(ver)[:]
q = np.ones(nboom - 2 + nspar)
qs1list = []
qs2list = []
qs3list = []

for pos in range(np.size(sc)):  
    xtop, xbot, ytop, ybot, Bup, Blow, xna, yna = \
    WB.Geometry(tskin, Aflange, ver, pos, sc)[0:8]
    
    qb, qbuplist, qblowlist = WB.openShearflow(sx, sy, ver, tskin, Aflange, pos, sc)[:]   
    
    L, ctbtl, ctbtm, ctbtr, eq00, eq01, eqr0, qs3to2, excess = \
        WB.rateoftwist(G, t, ver, tskin, Aflange, pos, sc, sx, sy)[:]
        
    qs1list, qs2list, qs3list, q, tau = \
    WB.qtot(sx, sy, pos, G, t, ver, tskin, Aflange, sc, q, qs1list, qs2list, qs3list)[:] 
    
    dthetadz = \
    WB.dthetadz(pos, G, t, ver, tskin, Aflange, sc, sx, sy, q, qs1list, qs2list, qs3list)
print q



"""tests"""
if ver == '6/1' or ver == '8/1':
    """up is positive"""
    Fy = qb[2*pos_spar1 - 1]*L[2*pos_spar1 - 1] - qb[pos_spar1 - 1]*L[pos_spar1 - 1] - \
            qb[nboom/2 + pos_spar1]*L[nboom/2 + pos_spar1]
elif ver == '6/2' or ver == '8/2':
    Fy = qb[2*pos_spar1 - 1]*L[2*pos_spar1 - 1] - qb[pos_spar1 - 1]*L[pos_spar1 - 1] - \
            qb[3*pos_spar1]*L[3*pos_spar1] - qb[nboom - 1]*L[nboom - 1]
            
#init number of Fx
#Fx = 0.
#for i in range(nboom - 2 + nspar):
#    """right is positive"""
#    if ver == '6/1' or ver == '8/1':
#        if i != pos_spar1 - 1 or i != pos_spar1*2 - 1 or i != nboom/2 + pos_spar1:
#            Fx = Fx + qb[i]*L[i]
#    elif ver == '6/2' or ver == '8/2':
#        if i != pos_spar1 - 1 or i != pos_spar1*2 - 1 or i != pos_spar1*3 or \
#        i != nboom/2 + pos_spar1:
#            Fx = Fx + qb[i]*L[i]
            
print 'Fy', Fy, 'Sy', sy
#print 'Fx', Fx, 'Sx', sx

#"""unit test"""
#for i in range(nboom + 2):
    
    
    
#    xtop, xbot, ytop, ybot, Bup, Blow, xna, yna, bup, blow, Al, Am, Ar = \
#    WB.Geometry(tskin, Aflange, ver, pos, sc)[:]
#        
#    Ixx, Iyy, Ixy = WB.MOI(tskin, Aflange, ver, pos, sc)[:]
#    
#    sigz, x, y, sigzbot, x_dlower, y_dlower = \
#    WB.NormalStress(i, sc, pos, tskin, Aflange, ver, Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[:]
#    
#    qb = WB.openShearflow(sx, sy, ver, tskin, Aflange, pos, sc)
#    
#    L, ctbtl, ctbtm, ctbtr, eq00, eq01, eqr0, qs3to2, excess = \
#    WB.rateoftwist(G, t, ver, tskin, Aflange, pos, sc, sx, sy)[:]    
#    
#    eq10, eq11, eqr1 = WB.torsionequivalence(sx, sy, pos, G, t, ver, tskin, Aflange, sc)[:]
    
#"""normal stress for the 2d plot"""
#    for j in range(8):
#        #stress on top panel        
#        z[pos, j] = NormalStress(pos, start + j*pitch_coor, start, pitch_coor, sc)[0]
#        #stress on bottom panel
#        zbot[pos, j] = NormalStress(pos, start + j*pitch_coor, start, pitch_coor, sc)[3]
        
#    """normal stress for 3d plot"""
#    for i in range(np.size(xup)):
#        #upper part
#        x = NormalStress(pos, i, start, pitch_coor, sc)[1]
#        y = NormalStress(pos, i, start, pitch_coor, sc)[2]
#        xylist[0].append(x*10.)
#        xylist[1].append(y*10.)
#        xylist[2].append(pos)
#        sigz = NormalStress(pos, i, start, pitch_coor, sc)[0]
#        xylist[3].append(sigz)
#        
#        #add lower part together
#        x_dlower = NormalStress(pos, i, start, pitch_coor, sc)[4]
#        y_dlower = NormalStress(pos, i, start, pitch_coor, sc)[-1]
#        xylist[0].append(x_dlower*10.)
#        xylist[1].append(y_dlower*10.)
#        xylist[2].append(pos)
#        sigzbot = NormalStress(pos, i, start, pitch_coor, sc)[3]
#        xylist[3].append(sigzbot)
        
#    """shear flow calc of each section"""
#    #get coefficients to solve qs1, qs1 using linalg.solve
#    
#    m = np.array([[eq00, eq01],[eq10, eq11]])
#    n = np.array([eqr0, eqr1])
#    qs1, qs2 = np.linalg.solve(m, n)[:]
#    """this calc qs1 and qs2, qs3 will be determined by qs2*qs3to2 + excess"""
#    qs3 = qs2*qs3to2 + excess
#    qs2list.append(qs2)
#    print qs1, qs2, qs3
#    #print qb[:2]
    
    
 
#    q0301, q0116, q1614 = \
#    openShearflow(sx, sy, pitch_coor, pos, pos_spar, nspar)[:3] #+ qs1
#    
#    qlist[0].append(q0301)
#    qlist[1].append(q0116)
#    qlist[2].append(q1614)    
#    
#    q1403 = openShearflow(sx, sy, pitch_coor, pos, pos_spar, nspar)[3] #+ qs1 - qs2
#    
#    qlist[3].append(q1403)
#    
#    q1409, q0908, q0803 = \
#    openShearflow(sx, sy, pitch_coor, pos, pos_spar, nspar)[4:] #+ qs2
#    
#    qlist[4].append(q1409)
#    qlist[5].append(q0908)
#    qlist[6].append(q0803)        
        
#"""Plot 2d bending stress"""
###color map plot
##fig = plt.figure()
##
###define x(i), and y(pos)
##pos = np.linspace(0, np.size(sc) - 2, np.size(sc) - 1)
##i = np.linspace(0, 7, 8)
##
###set colorbars to the same
##vmax1 = np.max(z)
##vmax2 = np.max(zbot)
##vmin1 = np.min(z)
##vmin2 = np.min(zbot)
##
##norm = mpl.colors.Normalize(min(vmin1, vmin2), max(vmax1, vmax2))
##
##sigmap = fig.add_subplot(121)
###gotta flip i and pos for some reason
##sigmap = pylab.pcolor(i, pos, z, norm=norm)
##
##pylab.colorbar(sigmap)
##sigmap.colorbar.set_label('Normal stress of top panel [Pa]')
##
##pylab.xlabel('Boom sections')
##pylab.ylabel('Span sections root to tip')
##
##sigbotmap = fig.add_subplot(122)
##sigbotmap = pylab.pcolor(i, pos, zbot, norm = norm)
##
##pylab.colorbar(sigbotmap)
##sigbotmap.colorbar.set_label('Normal stress fo bottom panel [Pa]')
##
##pylab.xlabel('Boom sections')
##pylab.ylabel('Span sections root to tip')
#
#"""3d plot of wing"""
##def plot_mayavi(data):
##    fig = mlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0))
##    obj = mlab.points3d(data[0], data[1], data[2], data[3], line_width = 1.0, mode="cube", scale_mode="none",scale_factor=0.5)
##    bar = mlab.colorbar(object=obj, title="Stress (in Pa)", orientation="vertical")
##    mlab.view(50,50)
##    mlab.roll(5)
##    #mlab.savefig(os.path.dirname(os.path.realpath(__file__))+'\\sim_outputs\\plot_mayavi.jpg', size=[5000,3000])
##    mlab.show()  
##    
##plot_mayavi(xylist)   
##pylab.show()
#
#"""plot wingbox geometry"""
#"""neutral axis not correctly dispalced"""
##plt.rcParams.update({'font.size': 16})
##
###preset pos for plotting purpose
##pos = 0
##sc[pos]=1.
###distance from boom to y = 0
##ytop = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[2]
###this is negative already
##ybot = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[3]
##
###neutral axis
##xna = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[6]
##yna = Geometry(nboom, tskin, pitch_coor, Aflange, pos, start, stop, pos_spar)[7]
##
###init x1 list, which is x axis
##x1 = np.zeros(nboom/2)
##for i in range(nboom/2):
##    x1[i] = x0[start + i*pitch_coor]*sc[pos]
##
##plt.plot(x1, ytop, linewidth = 3, linestyle = '-', c = 'r', label = 'Wingbox')
##plt.plot(x1, -ybot, linewidth = 3, linestyle = '-', c = 'r')
##
###plot neutral axis
##plt.plot(np.ones(nboom/2)*xna, np.linspace(-0.2, 0.2, nboom/2), 'b--', label = 'Neutral Axis')
##plt.plot(np.linspace(-0.1, 1.1, nboom/2), np.ones(nboom/2)*yna, 'b--')
##
###plot start and end spars
##plt.plot(np.ones(nboom/2)*x1[0], np.linspace(-ybot[0], ytop[0], nboom/2), \
##        linewidth = 3, linestyle = '-', c = 'r')
##plt.plot(np.ones(nboom/2)*x1[-1], np.linspace(-ybot[-1], ytop[-1], nboom/2), \
##        linewidth = 3, linestyle = '-', c = 'r')
##        
###plot middle spar which is at 4th boom from left
##plt.plot(np.ones(nboom/2)*x1[2], np.linspace(-ybot[2], ytop[2], nboom/2), \
##        linewidth = 3, linestyle = '-', c = 'r')
##
###plot airfoil
##plt.plot(sc[pos]*xup, sc[pos]*yup, 'g-')
##plt.plot(sc[pos]*xlow, sc[pos]*ylow, 'g-')
##plt.axis('equal')
##
###plot booms
##i = 0
##while i < nboom/2:     
##    circletop=plt.Circle((x1[i], ytop[i]), .01,color='r', clip_on = False)
##    fig = plt.gcf()
##    fig.gca().add_artist(circletop)
##    circletop=plt.Circle((x1[i], -ybot[i]), .01,color='r', clip_on = False)
##    fig = plt.gcf()
##    fig.gca().add_artist(circletop)
##    i = i + 1
##
##plt.xlabel('Chord length [m]')
##plt.ylabel('Camber [m]')
##plt.legend(loc = 1)
#
#"""plot shear flow diagram"""
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
#Z2 = np.ones(10)*qlist[2][0]
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
#"""dotted lines to connect shadow and q magnitude"""
##01 to q0116
#X010 = np.ones(10)*xup[start]
#Y010 = np.ones(10)*yup[start]
#Z010 = np.linspace(0, qlist[1][0], 10)
#ax.plot_wireframe(X010, Y010, Z010, color = 'g', linestyles = ':')
#
##16 to q0116
#X160 = np.ones(10)*xlow[start]
#Y160 = np.ones(10)*ylow[start]
#Z160 = np.linspace(0, qlist[1][0], 10)
#ax.plot_wireframe(X160, Y160, Z160, color = 'g', linestyles = ':')
#
##14 to q1614
#X140 = np.ones(10)*xlow[start + pitch_coor*(pos_spar - 1)]
#Y140 = np.ones(10)*ylow[start + pitch_coor*(pos_spar - 1)]
#Z140 = np.linspace(0, -qlist[2][0], 10)
#ax.plot_wireframe(X140, Y140, Z140, color = 'g', linestyles = ':')
#
##14 to q1409
#X141 = np.ones(10)*xlow[start + pitch_coor*(pos_spar - 1)]
#Y141 = np.ones(10)*ylow[start + pitch_coor*(pos_spar - 1)]
#Z141 = np.linspace(0, qlist[4][0], 10)
#ax.plot_wireframe(X141, Y141, Z141, color = 'r', linestyles = ':')
#
##03 to q1403
#X030 = np.ones(10)*xup[start + pitch_coor*(pos_spar - 1)]
#Y030 = np.ones(10)*yup[start + pitch_coor*(pos_spar - 1)]
#Z030 = np.linspace(0, qlist[3][0], 10)
#ax.plot_wireframe(X030, Y030, Z030, color = 'r', linestyles = ':')
#
##03 to q0301
#X031 = np.ones(10)*xup[start + pitch_coor*(pos_spar - 1)]
#Y031 = np.ones(10)*yup[start + pitch_coor*(pos_spar - 1)]
#Z031 = np.linspace(0, qlist[0][0], 10)
#ax.plot_wireframe(X031, Y031, Z031, color = 'g', linestyles = ':')
#
##09 to q1409
#X090 = np.ones(10)*xlow[stop]
#Y090 = np.ones(10)*ylow[stop]
#Z090 = np.linspace(0, qlist[4][0], 10)
#ax.plot_wireframe(X090, Y090, Z090, color = 'r', linestyles = ':')
#
##08 to q0908
#X080 = np.ones(10)*xup[stop]
#Y080 = np.ones(10)*yup[stop]
#Z080 = np.linspace(0, qlist[5][0], 10)
#ax.plot_wireframe(X080, Y080, Z080, color = 'r', linestyles = ':')
#
#ax.xaxis.set_major_formatter(plt.NullFormatter())
#ax.yaxis.set_major_formatter(plt.NullFormatter())
#ax.set_xlabel('Chord location [-]')
#ax.set_ylabel('Camber location [-]')
#ax.set_zlabel('Shear flow [N/m]')
#
##plt.tight_layout()
#plt.show()