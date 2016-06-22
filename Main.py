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
import Operation_Loading as w
from Configuration import *
#I needed to make __init__.py in Launch folder to make this work
import Launch.LaunchLoad as LL

loading = raw_input('Choose loading case: a.Launch; b.Cruise; c.Ascent; d.Descent Climbout; \
e.Entry; f.Turning: ')

#material property
"""PETI-330 (59%)"""
Sig_ten = 446.4*10**6/1.5
Sig_comp = 312.*10**6/1.5
E = 64.*10**9
Tau = 52.*10**6/1.5
"""not correct yet!"""
Pratio = 0.1

#from wing.py
b=50.
Cl=0.5796
Cd=0.024
rho = 1300.
V=40.       #velocity in m/s
dynaP=1000.
A=10.
taper=0.55
F=0.
tw=0.005
g=8.8
S=35.
rhowing=1300.
span = 8.106*2

taper1 = 0.55
taper2 = 0.55
S1 = 35.
S2 = 35.

"""need to change this to cfrp"""
G = 5.*10**9

#choose launcher type
launch_type = 'block 1'
#real payload mass
W_sc = 2000. #kg
#weight of a wing
"""should get this from WingboxBoom"""
W_wing = 185. #kg

hingel = 1.4#span/2.

"""from LaunchLoad"""
L_fs = 5.
W_pl = 1000.
R = 0.6
t_fs = 0.005

#number of sections for analysis
n = 200
#launch
if loading == 'a':
    """this is not correct atm"""
#    #sx is drag direction force
#    sx = LL.WingLaunchStress(W_wing, launch_type, W_sc)[4]#/float(n)*np.ones(n)
#    #sy is lift direction force
#    sy = LL.WingLaunchStress(W_wing, launch_type, W_sc)[5]#/float(n)*np.ones(n)
    
    sx = w.loadcase(Cl, Cd, dynaP, A, taper1, taper2, S1, S2, rhowing, F, tw, \
    g, n, loading, hingel, span)[1]
    sy = w.loadcase(Cl, Cd, dynaP, A, taper1, taper2, S1, S2, rhowing, F, tw, \
    g, n, loading, hingel, span)[0]

#cruise
elif loading  == 'b' or loading == 'f':
    Cl = 0.5796
    Cd = 0.024
    sx = w.loadcase(Cl, Cd, dynaP, A, taper1, taper2, S1, S2, rhowing, F, tw, \
    g, n, loading, hingel, span)[1]
    sy = w.loadcase(Cl, Cd, dynaP, A, taper1, taper2, S1, S2, rhowing, F, tw, \
    g, n, loading, hingel, span)[0]

#ascent
elif loading == 'c':
    """need confirmation"""
    Cl = 0.18
    Cd = 0.008
    sx = w.loadcase(Cl, Cd, dynaP, A, taper1, taper2, S1, S2, rhowing, F, tw, \
    g, n, loading, hingel, span)[1]
    sy = w.loadcase(Cl, Cd, dynaP, A, taper1, taper2, S1, S2, rhowing, F, tw, \
    g, n, loading, hingel, span)[0] 
    

#descent
elif loading == 'd':
    Cl = 0.06
    Cd = 0.005
    sx = w.loadcase(Cl, Cd, dynaP, A, taper1, taper2, S1, S2, rhowing, F, tw, \
    g, n, loading, hingel, span)[1]
    sy = w.loadcase(Cl, Cd, dynaP, A, taper1, taper2, S1, S2, rhowing, F, tw, \
    g, n, loading, hingel, span)[0]
    

#entry
elif loading == 'e':
    Cl = 1.
    Cd = 1.
    sx = w.loadcase(Cl, Cd, dynaP, A, taper1, taper2, S1, S2, rhowing, F, tw, \
    g, n, loading, hingel, span)[1]
    sy = w.loadcase(Cl, Cd, dynaP, A, taper1, taper2, S1, S2, rhowing, F, tw, \
    g, n, loading, hingel, span)[0] 
    
    
else:
    print 'Really? Do it again!'

"""Parameters. Change this if needed!"""
        
#tskin = 0.005
#Aflange = 300*10**(-6)




#scaling of airfoil
sc = w.loadcase(Cl, Cd, dynaP, A, taper1, taper2, S1, S2, rhowing, F, tw, g, \
n, loading, hingel, span)[-1]

sectionlist = ['roottohinge', 'hingetosec2', 'sec2totip']
verlist = ['4/1', '4/2', '5/1', '5/2', '6/1', '6/2', '7/1', '7/2', '8/1', '8/2']

    

#pos = 5
#Mx = (w.loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[2])[pos]
#print Mx

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


"""commented this out cuz conflict with normal stress"""
#for pos in range(np.size(sc)):  
#    xtop, xbot, ytop, ybot, Bup, Blow, xna, yna = \
#    WB.Geometry(tskin, Aflange, ver, pos, sc)[0:8]
#    
#    qb, qbuplist, qblowlist = WB.openShearflow(sx, sy, ver, tskin, Aflange, pos, sc)[:]   
#    
#    L, ctbtl, ctbtm, ctbtr, eq00, eq01, eqr0, qs3to2, excess = \
#    WB.rateoftwist(G, t, ver, tskin, Aflange, pos, sc, sx, sy)[:]
        
#    qs1list, qs2list, qs3list, q, tau = \
#    WB.qtot(sx, sy, pos, G, t, ver, tskin, Aflange, sc, q, qs1list, qs2list, qs3list)[:] 
#    
#    dthetadz = \
#    WB.dthetadz(pos, G, t, ver, tskin, Aflange, sc, sx, sy, q, qs1list, qs2list, qs3list)

#init dimension array full of inf, every run a row will be filled with tskin, section, Aflange, ver
dimension = np.full((0, 7), np.inf)

#init mass empty list
#append add a list/number directly to the back of existing list (same row)
#extend add elements of a list to the back of existing list (same row)
masslist = []
for i in range(len(verlist)):
    masslist.append([])

tskinlist = np.arange(0.0005, 0.003, 0.0015)
tskinlist = tskinlist.tolist()

tsparlist = np.arange(0.0005, 0.003, 0.0015)
tsparlist = tsparlist.tolist()

Aflangelist = np.arange(50.*10**(-6), 350.*10**(-6), 150*10**(-6))
Aflangelist = Aflangelist.tolist()

#for tskin in tskinlist:
#    #to print how many percent finished
#    print '' 
#    print '---', float(tskinlist.index(tskin))/len(tskinlist)*100, '% fninished---'
#    print ''
#    
#    for tspar in tsparlist:        
#        #section after thickness because I don't want different t in different sections
#        for section in sectionlist:
#        #section = 'roottohinge'
#        
#            for Aflange in Aflangelist:
#    #        #get configuration parameters
#    #        start, stop, nboom, pitch, pitch_coor, pos_spar1, pos_spar2, nspar = \
#    #        config(ver)
#            
#            #make a section selection loop
#                for ver in verlist:
#                #ver = '4/1'
#                    #do this so that large hopeless values are discarded
#                    if ver == '8/1' or ver == '8/2' or ver == '7/1' or ver == '7/2':
#                        if float(tskinlist.index(tskin))/len(tskinlist) > 0.5 and \
#                        float(tsparlist.index(tspar))/len(tsparlist) > 0.5 and \
#                        float(Aflangelist.index(Aflange))/len(Aflangelist) > 0.5:
#                            continue
#                
#                    print 'Section:', section, 'Version:', ver, 'tskin', tskin, \
#                    'tspar', tspar, 'Aflange', Aflange
#                    
#                
#                    #get config from config def
#                    start, stop, nboom, pitch, pitch_coor, pitch_coor_last, \
#                    pos_spar1, pos_spar2, nspar = config(ver, section)[:9]
#                    
#                    q = np.ones(nboom - 2 + nspar)
#                    
#                    pos_t, pos_r, L_sec = spanConfig(section, span, sc, hingel)
#                    """The estimation starts here for each section"""
#                    """all lists refresh here"""
#                    
#        
#                    qs1list = []
#                    qs2list = []
#                    qs3list = [] 
#                    
#                    Tsaihill_toplist = []
#                    Tsaihill_botlist = []
#                    
#                    #to get the avg area of ribs, avg area of booms in a section
#                    pos_avg = pos_t/2
#                    
#                    #get spar height, boom area for weight calc
#                    ytop, ybot, Bup, Blow = WB.Geometry(tskin, Aflange, ver, \
#                    pos_avg, sc, section, tspar)[2:6]
#                    #get boom distance, cell area for weight calc at each position
#                    bup, blow, Al, Am, Ar = WB.Geometry(tskin, Aflange, ver, \
#                    pos_avg, sc, section, tspar)[-5:]   
#                    
#    #                if section == 'sec2totip' or section == 'hingetosec2':
#    #                    if ver == '8/2' or ver == '7/2' or ver == '6/2' or ver == \
#    #                    '5/2' or ver == '4/2':
#    #                        sigz  = 50*26.1*9.81/(sum(Bup) + sum(Blow) + \
#    #                        tspar*(ytop[3] + ybot[3])*4)
#    #                    else:
#    #                        sigz  = 50*26.1*9.81/(sum(Bup) + sum(Blow) + \
#    #                        tspar*(ytop[3] + ybot[3])*3)
#                    
#                    for pos in np.arange(0, 11, 1):
#                        #refresh normal stress list so every time siglist caps stress 
#                        #across a certain chord
#                        sigzlist = []
#                        sigzbotlist = []
#                        
#                        for i in range(151):
#                #                Mx = (w.loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[2])[pos]
#                #                My = (w.loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[3])[pos]
#                        
#                            sigz, x, y, sigzbot, x_dlower, y_dlower = \
#                            WB.NormalStress\
#                            (i, sc, tskin, Aflange, ver, Cl, Cd, dynaP, A, taper1, \
#                            taper2, S1, S2, rhowing, F, tw, g, span, section, hingel, \
#                            n, loading, tspar, pos)
#                            
#                            sigzlist.append(sigz)
#                            sigzbotlist.append(sigzbot)
#                            
#                        sigzmax = max(np.abs(sigzlist))
#                        sigzbotmax = max(np.abs(sigzbotlist))
#                               
#                        
#                        
#                        qs1list, qs2list, qs3list, q, tau = \
#                        WB.qtot(sx, sy, pos, G, ver, tskin, Aflange, sc, q, qs1list, \
#                        qs2list, qs3list, section, tspar)[:]
#                        
#                        #use np.abs here to get abs of a list, which can't be done with abs
#                        taumax = max(np.abs(tau))
#                        
#                        #Tsai Hill value for top panel
#                        Tsaihill_top = sigzmax**2/Sig_ten**2 - sigzmax*sigzmax*\
#                        Pratio/Sig_ten**2 + (sigzmax*Pratio)**2/Sig_comp**2 + \
#                        taumax**2/Tau**2
#                        
#                        Tsaihill_toplist.append(Tsaihill_top)
#                        
#                        #Tsai Hill value for bot panel
#                        Tsaihill_bot = sigzbotmax**2/Sig_ten**2 - sigzbotmax*\
#                        sigzbotmax*Pratio/Sig_ten**2 + (sigzbotmax*Pratio)**2/\
#                        Sig_comp**2 + taumax**2/Tau**2
#                        
#                        Tsaihill_botlist.append(Tsaihill_bot)
#                        
#    #                Tsaihill_toplist = [0]
#    #                Tsaihill_botlist = [0]
#    
#                    print 'Checking'
#                    print ''
#                    
#                    #Tsai Hill failure creteria
#                    if (max(Tsaihill_toplist) < 1.) and (max(Tsaihill_botlist) < 1.):
#                        
#                        print 'Section:', section, 'Version:', ver
#                        print 'tskin, tspar, Aflange:', tskin, tspar, Aflange
#                        jlist = WB.buckling(span, tskin, E, i, sc, Aflange, ver, \
#                        sigzmax, section, hingel, tspar)
#                        
#                        #set up a filter to avoid the case when all spar number are unmatched
#                        if jlist != []:                    
#                        
#                            
#                            if ver == '8/2' or ver == '7/2' or ver == '6/2' or \
#                            ver == '5/2' or ver == '4/2':
#                                #calc mass with least number of spars
#                                mass = ((sum(Bup) + sum(Blow))*L_sec + \
#                                tskin*(Al + Am + Ar)*min(jlist) + (sum(bup) + sum(blow))*\
#                                L_sec*tskin + tspar*(ytop[3] + ybot[3])*4)*rho
#                            else:
#                                mass = ((sum(Bup) + sum(Blow))*L_sec + \
#                                tskin*(Al + Am + Ar)*min(jlist) + (sum(bup) + sum(blow))*\
#                                L_sec*tskin + tspar*(ytop[3] + ybot[3])*3)*rho
#                            
#                            #add corresponding parameters to a 2d array
#                            dimension = np.vstack((dimension, [ver, section, tskin, \
#                            tspar, Aflange, min(jlist), mass]))
#                            
#                            #add elements to a list for mass calc
#                            masslist[verlist.index(ver)].extend([mass])
#                            
#                            print ''
#                            
#                            #print 'dimension', dimension
#                            print 'Mass list expanding:'
#                            print masslist
#                            print ''
                        
                        

#not necessary anymore since dimension is an empty array to begin with                    
#dimension = dimension[1:][:]

'''potentially use this from nasa: Efficiency, and Design Data for Beryllium
Structures.'''
#e = 0.556
#Wbend = sy.Symbol('Wbend')
#Wbend = sy.solve(Wbend/rho/Zs/t - eps*(M/Zs/t**2/E)**e, Wbend)
#Wshear = rho*Fs/tau*1.5

"""for moment diagram only"""
"""cruise"""
#moment around x and moment around y
#Mx = (w.loadcase(Cl, Cd, dynaP, A, taper1, taper2, S1, S2, rhowing, F, tw, \
#g, n, loading, hingel, span)[2])
#My = (w.loadcase(Cl, Cd, dynaP, A, taper1, taper2, S1, S2, rhowing, F, tw, \
#g, n, loading, hingel, span)[3])
#
#fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()
#plt.rcParams.update({'font.size': 16})
#
#axis1 = ax1.plot(sc*3.178823529411765 - 4.609294117647059, sy, 'g', label = '$V_y$')
#axis2 = ax1.plot(sc*3.178823529411765 - 4.609294117647059, sx, 'g--', label = '$V_x$')
#
#if loading == 'e':
#    vz = np.zeros(36)
#    vz = vz.tolist()
#    for i in range(len(sc) - 36):
#        
#        vz.extend([185./4*26.1*9.81*np.cos(0.1611)])
#    axis3 = ax1.plot(sc*3.178823529411765 - 4.609294117647059, vz, 'g-.', label = '$V_z$')
#    
#elif loading == 'a':
#    vz = np.ones(36)*62.8923326411827
#    vz = vz.tolist()
#    for i in range(len(sc) - 36):
#        
#        vz.extend([5995.87683575586/4*np.cos(0.1611)])
#    axis3 = ax1.plot(sc*3.178823529411765 - 4.609294117647059, vz, 'g-.', label = '$V_z$')
#
##axis2 = ax2.plot(sc*3.178823529411765 - 4.609294117647059, Mx, 'b--', label = '$M_x$')
#axis4 = ax2.plot(sc*3.178823529411765 - 4.609294117647059, My, 'b:', label = '$M_x = M_y = 0$')
##ax1.set_ylim(ymin = -850)
#ax1.set_xlim(xmin = 0)
##ax2.set_ylim(ymin = 0)
#ax1.set_xlabel('Wing span from tip to root [$m$]')
#
##vertical axis at x = 8.106
#plt.axvline(8.106, color = 'r', linestyle = ':')
#
#ax1.set_ylabel(r'Loading [$N$]', color = 'g')
#ax2.set_ylabel(r'Moment [$Nm$]', color = 'b')
#
#lns = axis1 + axis2 +axis3 + axis4
#labs = [l.get_label() for l in lns]
#ax1.legend(lns, labs, loc='upper left', shadow = True)
#
##ax1.legend(loc="upper left", shadow=True, fancybox=True)
##           
##ax2.legend(loc="lower left", shadow=True, fancybox=True)
#
#ax1.annotate('Root line', xy=(8.106, 3000), xytext=(7.2, 4000),
#            arrowprops=dict(arrowstyle="fancy", #linestyle="dashed",
#                            color="0.5",
#                            #patchB=el,
#                            shrinkB=5,
#                            connectionstyle="arc3,rad=0.3",
#                            ),
#                        )
#
#plt.tight_layout()
#plt.show()



"""tests"""
#if ver == '6/1' or ver == '8/1':
#    """up is positive"""
#    Fy = qb[2*pos_spar1 - 1]*L[2*pos_spar1 - 1] - qb[pos_spar1 - 1]*L[pos_spar1 - 1] - \
#            qb[nboom/2 + pos_spar1]*L[nboom/2 + pos_spar1]
#elif ver == '6/2' or ver == '8/2':
#    Fy = qb[2*pos_spar1 - 1]*L[2*pos_spar1 - 1] - qb[pos_spar1 - 1]*L[pos_spar1 - 1] - \
#            qb[3*pos_spar1]*L[3*pos_spar1] - qb[nboom - 1]*L[nboom - 1]
#            
##init number of Fx
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
            
#print 'Fy', Fy, 'Sy', sy
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
"""plot wingbox geometry"""
"""neutral axis not correctly dispalced"""
plt.rcParams.update({'font.size': 16})

#preset pos for plotting purpose
pos = 0
sc[pos]=1.

tskin = 0.0005
Aflange = 100*10**(-6)
tspar = 0.001

for ver in verlist[:4]:
    #ver = '8/2'
    section = 'roottohinge'
    
    
    start, stop, nboom, pitch, pitch_coor, pitch_coor_last, pos_spar1, \
                pos_spar2, nspar= config(ver, section)[:9]
    
    #coor of booms
    xtop, xbot, ytop, ybot = WB.Geometry(tskin, Aflange, ver, pos, sc, \
    section, tspar)[:4]
    
    #neutral axis
    xna = WB.Geometry(tskin, Aflange, ver, pos, sc, section, tspar)[6]
    yna = WB.Geometry(tskin, Aflange, ver, pos, sc, section, tspar)[7]
    
#    #init x1 list, which is x axis
#    x1 = np.zeros(nboom/2)
#    for i in range(nboom/2):
#        x1[i] = xraw[start + i*pitch_coor]*sc[pos]
        
    plt.subplot(2, 2, verlist.index(ver) + 1)
    plt.plot(xtop, ytop, linewidth = 3, linestyle = '-', c = 'r', label = 'Wingbox')
    plt.plot(xtop, -ybot, linewidth = 3, linestyle = '-', c = 'r')
    
    #plot neutral axis
    plt.plot(np.ones(nboom/2)*xna, np.linspace(-0.2, 0.2, nboom/2), 'b--', label = 'Neutral Axis')
    plt.plot(np.linspace(-0.1, 1.1, nboom/2), np.ones(nboom/2)*yna, 'b--')
    
    #plot start and end spars
    plt.plot(np.ones(nboom/2)*xtop[0], np.linspace(-ybot[0], ytop[0], nboom/2), \
            linewidth = 3, linestyle = '-', c = 'r')
    plt.plot(np.ones(nboom/2)*xtop[-1], np.linspace(-ybot[-1], ytop[-1], nboom/2), \
            linewidth = 3, linestyle = '-', c = 'r')
            
    #plot 1st spar which is at 3rd/2nd boom from left
    plt.plot(np.ones(nboom/2)*xtop[pos_spar1 - 1], np.linspace(-ybot[pos_spar1 - 1], \
    ytop[pos_spar1 - 1], nboom/2), \
            linewidth = 3, linestyle = '-', c = 'r')
            
    """plot 2nd spar which is at pos_spar2 - 1"""
    if ver == '4/2' or ver == '5/2' or ver == '6/2' or ver == '7/2' or ver == '8/2':
        plt.plot(np.ones(nboom/2)*xtop[pos_spar2 - 1], np.linspace(-ybot[pos_spar2 - 1], \
    ytop[pos_spar2 - 1], nboom/2), \
            linewidth = 3, linestyle = '-', c = 'r')
    
    #plot airfoil
    plt.plot(sc[pos]*xup, sc[pos]*yup, 'g-')
    plt.plot(sc[pos]*xlow, sc[pos]*ylow, 'g-')
    plt.axis('equal')
    
    #plot booms
    i = 0
    while i < nboom/2:     
        circletop=plt.Circle((xtop[i], ytop[i]), .01,color='r', clip_on = False)
        fig = plt.gcf()
        fig.gca().add_artist(circletop)
        circletop=plt.Circle((xbot[i], -ybot[i]), .01,color='r', clip_on = False)
        fig = plt.gcf()
        fig.gca().add_artist(circletop)
        i = i + 1
    
    plt.xlabel('Chord length [m]')
    plt.ylabel('Camber [m]')
    
    """x, y limit doesn't due to plt.axis('equal')"""
    axes = plt.gca()
    axes.set_xlim([0,1])
    #axes.set_ylim([-0.3, 0.3])
    
    plt.legend(loc = 1)

plt.tight_layout()
plt.show()
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
