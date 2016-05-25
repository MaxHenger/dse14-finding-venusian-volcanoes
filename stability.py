# -*- coding: utf-8 -*-
"""
Created on Wed May 18 15:34:27 2016

@author: Julius
"""
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

#from pylab import *
from matplotlib.widgets import Slider
#import wing

class __dummy_wing__:
    def __init__(self):
        self.surface=1.
        self.x_location=1.
        self.dist_np=1.
        self.wash=0.
        self.VelFrac=1.
        self.clalpha=1.
        self.chord=1.
        self.cm = 0.
        self.cl= 1.
    def __str__(self):
        print(self.surface)
        print(self.dist_np)
        print(self.wash)
        print(self.VelFrac)
        print(self.clalpha)
        print(self.chord)

canard = __dummy_wing__()
canard.surface=2.
canard.dist_np=-7.5
canard.clalpha=0.0705
canard.wash=0
canard.VelFrac=1
canard.cl=0.23

main = __dummy_wing__()
main.surface=35.3
main.cl = 0.23
main.cm = -0.093
main.chord = 2.633
main.clalpha = 0.0705

tail = __dummy_wing__()
tail.surface=0
tail.dist_np=7.5
tail.cl = -0.23
tail.clalpha=0.0705
tail.wash = 0.1
tail.VelFrac = 0.8


xcg = 2.
xac = 1.23
stabMargin=0.1

# with cg and np input return sizing of tail / canard
def sizeStab(xcg,xac,canard,main,tail,stabMargin=0.1,ret="both"):
    """Returns area of either canard, tail, or coefficients of the equation ret= both, c, t
       Input: Canard, Wing, Tail location and aerodynamic properties
       Ouput: Canard, Tail surface Area"""
#    wing.surface # surface area
#    wing.x_location or wing.dist_np
#    wing.wash # upwash (positive) or downwash(negative) due to the main wing, for the main wing that can be 0 ( d eta / d cl_alpha_main )
#    wing.VelFrac # the velocity fraction with respect to the main wing
#    wing.clalpha # CL_alpha
#    wing.chord
#    wing.exist
    
    import sympy
    ScA=sympy.Symbol("ScA")
    StA=sympy.Symbol("StA")
    Sxcg=sympy.Symbol("Sxcg")
    
    eq = xac/main.chord - stabMargin - Sxcg/main.chord \
        + ScA*canard.dist_np/(main.surface*main.chord)*(canard.VelFrac)**2 *canard.clalpha/main.clalpha*(1+canard.wash)\
        + StA*tail.dist_np  /(main.surface*main.chord)*(tail.VelFrac)**2   *tail.clalpha  /main.clalpha*(1+tail.wash)
    
    if ret=="both":
        eqP = sympy.poly(eq)
        return eq,eqP.coeffs()
    elif ret=="c":
        return sympy.solve(eq.subs(StA,0),ScA)
    elif ret=="t":
        return sympy.solve(eq.subs(ScA,0),StA)
        
def sizeControl(xcg,xac,canard,main,tail,ret="both"):
    """Returns area of either canard, tail, or coefficients of the equation ret= both, c, t
       Input: Canard, Wing, Tail location and aerodynamic properties
       Ouput: Canard, Tail surface Area"""
    import sympy
    ScA=sympy.Symbol("ScA")
    StA=sympy.Symbol("StA")
    Sxcg=sympy.Symbol("Sxcg")
    
    # Sxcg => Symbolic X location of Center of Gravity
    # xac => Xlocation of aerodynamic center
    # 
    eq  = xac/main.chord - Sxcg/main.chord - main.cm/main.cl \
         + canard.cl/main.cl*ScA*canard.dist_np/(main.surface*main.chord)*(canard.VelFrac)**2 \
         -   tail.cl/main.cl*StA*  tail.dist_np/(main.surface*main.chord)*(  tail.VelFrac)**2
            
    #sympy.plotting.plot3d(Sxcg,(Sc,-5,5))
    if ret=="both":
        eqP = sympy.poly(eq)
        return eq,eqP.coeffs()
    elif ret=="c":
        return sympy.solve(eq.subs(StA,0),ScA)
    elif ret=="t":
        return sympy.solve(eq.subs(ScA,0),StA)
        
def plot_sizing():
    pass
#def manual_sizing(xcg,xac,canard,main,tail,stabMargin=0.05):
print("\n")
print(sizeStab(xcg,xac,canard,main,tail,stabMargin,ret="c"))
print(sizeControl(xcg,xac,canard,main,tail,ret="c"))




def variablePlot():
    variables=['canard_length','tail_length','surface','chord','Vcfrac','Vtfrac','washc','washt','cm_ac','CLc','CLt','CLah']
    init     =[canard.dist_np,tail.dist_np,main.surface,main.chord,canard.VelFrac,tail.VelFrac,canard.wash,tail.wash,main.cm,canard.cl,tail.cl,main.cl]
    value    = init
    lim_l    =[0            ,0     ,0     ,0    ,0     ,0     ,0    ,0, -1,0,-1,0.01]
    lim_u    =[10           ,10    ,100     ,10   ,1     ,1     ,1    ,1, 1,1,1,2]
    s_axes =[]
    slider   =[]
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    subplots_adjust(left=0.1, bottom=0.7)
    Sc = np.arange(0, 15, 1)
    St = np.arange(0, 15, 1)
    Sc, St = np.meshgrid(Sc, St)
    maxcg = xac - stabMargin \
    + Sc*canard.dist_np/(main.surface*main.chord)*(canard.VelFrac)**2 *canard.clalpha/main.clalpha*(1+canard.wash)\
    + St*tail.dist_np    /(main.surface*main.chord)*(tail.VelFrac)**2   *tail.clalpha  /main.clalpha*(1+tail.wash)

    maxcg2 = xac - value[8]/value[11] + value[9]/value[11]*Sc*value[0]/(value[2]*value[3])*(value[4])**2 \
            -value[10]/value[11]  * St*value[1]    /(value[2]*value[3])*(value[5])**2
    
    surf = ax.plot_surface(Sc, St, maxcg, rstride=1, cstride=1, color="red" ,linewidth=0, antialiased=False) #cmap=cm.coolwarm
    surf2 = ax.plot_surface(Sc, St, maxcg2, rstride=1, cstride=1, color="blue",linewidth=0, antialiased=False)
    
    axcolor = 'lightgoldenrodyellow'
    for i in range(0,len(variables)):
        s_axes.append( axes([0.25, 0.1+0.05*i, 0.65, 0.03], axisbg=axcolor))
        slider.append( Slider(s_axes[-1],variables[i],lim_l[i],lim_u[i],valinit=init[i]) )
    
    def update(val):
        ax.clear()
        Sc = np.arange(0, 15, 1)
        St = np.arange(0, 15, 1)
        Sc, St = np.meshgrid(Sc, St)
        maxcg = xac - stabMargin \
        + Sc*value[0]/(value[2]*value[3])*(value[4])**2 *canard.clalpha/main.clalpha*(1+value[6])\
        + St*value[1]    /(value[2]*value[3])*(value[5])**2   *tail.clalpha  /main.clalpha*(1+value[7])
        
        maxcg2 = xac - value[8]/value[11] + value[9]/value[11]*Sc*value[0]/(value[2]*value[3])*(value[4])**2 \
            -value[10]/value[11]  * St*value[1]    /(value[2]*value[3])*(value[5])**2        
        
        for i in range(0,len(variables)):
            value[i] = slider[i].val
        surf2 = ax.plot_surface(Sc, St, maxcg2, rstride=1, cstride=1, color="blue", alpha=.5,linewidth=0, antialiased=False) 
        surf = ax.plot_surface(Sc, St, maxcg, rstride=1, cstride=1, color="red",linewidth=0, antialiased=False)
        show()
        
    for i in slider:
        i.on_changed(update)

    show()
    
