# -*- coding: utf-8 -*-
"""
Created on Wed May 18 15:34:27 2016

@author: Julius
"""

import numpy as np
import sympy
import matplotlib.pyplot as plt
import AircraftStabilityCoeff
import Utility as util
import AircraftStability as ACStab
import control

class __dummy_wing__:
    def __init__(self):
        self.surface=1.
        self.x_location=1.
        self.dist_np=1.
        self.wash=0.
        self.aspect=None
        self.VelFrac=1.
        self.clalpha=1.
        self.taper=0.55
        self.root=0
        self.chord=0
        self.span=0
        self.oswald=0.9
        self.dist_z=0
        self.sweep=0
        self.cm = 0.
        self.cl= 1.
        self.cd=0
    def __str__(self):
        print(self.surface)
        print(self.dist_np)
        print(self.wash)
        print(self.VelFrac)
        print(self.clalpha)
        print(self.chord)

def dummyWings():
    canard = __dummy_wing__()
    #canard.surface=2.
    canard.dist_np=-3.75
    canard.clalpha=0.1#0.054480873504886229#0.0705
    canard.wash=0
    canard.VelFrac=1
    canard.cl=0.23
    canard.clde=-0.72192708
    canard.sweep=0 # deg
    canard.aspect=8
    canard.surface=7.2#8.8284814051829894
    
    main = __dummy_wing__()
    main.surface=36.5
    main.cl = 0.23
    main.cd = main.cl/7.
    main.cm = -0.093
    main.root = 4
    main.taper=0.4375  
    main.span = main.surface/(main.root/2*(1+main.taper))
    main.chord = main.root*2./3*(1+main.taper+main.taper**2)/(1+main.taper)
    main.sweep = np.rad2deg(np.arctan(main.root*(1-main.taper)/2./main.span)) # deg
    main.clalpha = 0.1#0.086793918127679101#0.0705
    main.dihedral=0*np.pi/180
    main.oswald=0.9
    main.aspect=main.surface/(main.root/2*(1+main.taper))**2
    
    tail = __dummy_wing__()
    #tail.surface=0
    tail.dist_np=5.5
    tail.cl = 0#-0.23
    tail.clde = 0.9
    tail.clalpha=0.1#0.054480873504886229#0.0705
    tail.wash = -0.043688739886377427
    tail.aspect=2.5
    tail.sweep=0 # deg
    tail.surface=14.656962810365979
    tail.VelFrac = 0.85 # .85 for fuselage mounted, .95 for fin mounted, 1.00 for canards and Tail
    
    vert = __dummy_wing__()
    vert.dist_np=tail.dist_np
    vert.chord  = 0
    vert.span   = 0
    vert.surface= 1.5*4.04#0.05*main.surface
    vert.aspect = 1.4#1.4 # 1.4
    vert.cl     = 0
    vert.clalpha= 0.1 # naca 0014
    vert.clde   = 0.9
    vert.VelFrac= 0.95
    vert.dist_z = 1
    vert.sweep  = 35*np.pi/180
    vert.wash   = -0.01
    return canard,main,tail,vert


def sizeStab(xac,canard,main,tail,stabMargin=0.1,configuration="both",ratio=0):
    """Returns area of either canard, tail, or coefficients of the equation ret= both, c, t
       Input: Canard, Wing, Tail location and aerodynamic properties
       Ouput: Canard, Tail surface Area"""
    
    ScA=sympy.Symbol("ScA")
    StA=sympy.Symbol("StA")
    Sxcg=sympy.Symbol("Sxcg")
    
    eq = xac/main.chord - stabMargin - Sxcg/main.chord \
        + ScA*canard.dist_np/(main.surface*main.chord)*(canard.VelFrac)**2 *canard.clalpha/main.clalpha*(1+canard.wash)\
        + StA*tail.dist_np  /(main.surface*main.chord)*(tail.VelFrac)**2   *tail.clalpha  /main.clalpha*(1+tail.wash)
    
    if configuration=="both":
        eqP = sympy.poly(eq)
        return eq,eqP.coeffs()
    elif configuration=="c":
        return sympy.solve(eq.subs({StA:ratio*ScA}))
    elif configuration=="t":
        return sympy.solve(eq.subs({ScA:ratio*StA}))
    else:
        raise ValueError("Configuration Unknown")
        
def sizeControl(xac,canard,main,tail,configuration="both",ratio=0):
    """Returns area of either canard, tail, or coefficients of the equation ret= both, c, t
       Input: Canard, Wing, Tail location and aerodynamic properties
       Ouput: Canard, Tail surface Area"""
    ScA=sympy.Symbol("ScA")
    StA=sympy.Symbol("StA")
    Sxcg=sympy.Symbol("Sxcg")
    
    eq  = (xac - Sxcg)/main.chord - main.cm/main.cl \
         + canard.cl/main.cl*ScA*canard.dist_np/(main.surface*main.chord)*(canard.VelFrac)**2 \
         +   tail.cl/main.cl*StA*  tail.dist_np/(main.surface*main.chord)*(  tail.VelFrac)**2
            
    if configuration=="both":
        eqP = sympy.poly(eq)
        return eq,eqP.coeffs()
    elif configuration=="c":
        return sympy.solve(eq.subs(ScA,ratio*StA))
    elif configuration=="t":
        return sympy.solve(eq.subs(StA,ratio*ScA))
    else:
        raise ValueError("Configuration Unknown")
        
    
def return_sizing(xac,canard,main,tail,configuration="t",ratio=0,safety=1,ran=(-2,10),step=0.001,xcg=0,plot=False):
    import sympy.plotting
    stab=sizeStab(xac,canard,main,tail,stabMargin,configuration,ratio)[0]
    cont=sizeControl(xac,canard,main,tail,configuration,ratio)[0]
    
    x=np.arange(ran[0],ran[1],step)
    f1 = sympy.lambdify(stab.get(stab.keys()[0]).free_symbols.pop(),stab.get(stab.keys()[0]))
    f2 = sympy.lambdify(cont.get(cont.keys()[0]).free_symbols.pop(),cont.get(cont.keys()[0]))
    y1=f1(x)
    y2=f2(x)
    S_l1 = f1(xcg[0])
    S_l2 = f2(xcg[0])
    S_r1 = f1(xcg[1])
    S_r2 = f2(xcg[1])
    S = max([S_l1,S_l2,S_r1,S_r2])
    
    if plot:
        fig = plt.figure("Stability and Control")
        ax=fig.add_subplot(111)
        ax.plot(x,y1,linestyle="-",linewidth=1,label=r"$\rm Stability$")
        ax.plot(x,y2,linestyle="--",linewidth=1.5,label=r"$\rm Control$")
        ax.axvline(xcg[0],color="r",linewidth=1.5,linestyle="-.",label=r"$\rm CG$")
        ax.axvline(xcg[1],color="r",linewidth=1.5,linestyle="-.")        
        ax.axhline(S,color="c",linewidth=1.,linestyle="--",label=r"$\rm S_{tail,min}$")
        ax.legend(loc=5)
        ax.grid(True)
        
        ax.set_title(r"Stability and Control", fontsize=14)
        ax.set_xlabel(r"x location", fontsize=14)
        ax.set_ylabel(r"Surface area tail", fontsize=14)
        ax.xaxis.set_label_coords(0.5, -0.025)
        ax.yaxis.set_label_coords(-0.025,0.5)
        #plt.figtext(0.30, 0.03, "configuration: "+configuration)
        #if configuration=="t":
        #    plt.figtext(0.30, 0.06, "S tail/canard ratio: "+str(ratio))
        #    ax.set_ylabel(r"Surface area tail", fontsize=14)
        #elif configuration=="c":
        #    plt.figtext(0.30, 0.06, "S canard/tail ratio: "+str(ratio))
        #    ax.set_ylabel(r"Surface area canard", fontsize=14)
        
        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position('zero')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_ylim([-5,50])
        fig.show()
    
    y=abs(y1-y2)
    x_min=x[y.argmin()]
    S_min = max([f1(x_min)*safety , f2(x_min)*safety])
    
    if configuration=="t":
        return x_min,S_min,S*ratio,S
    elif configuration=="c":
        return x_min,S_min,S,S*ratio
    
def cmalpha(canard,main,tail,xac,xcg):
    mw = 0 # dCm / d alpha
    mw = main.clalpha*(xcg-xac)/main.chord \
         - tail.clalpha*tail.dist_np*tail.surface/(main.surface*main.chord)*tail.VelFrac**2 \
         - canard.clalpha*canard.dist_np*canard.surface/(main.surface*main.chord)*canard.VelFrac**2
    
    return mw
    
def DATCOM_tail(mach,tail,etha=0.95):
    beta = np.sqrt(1-mach**2)
    CL = 2*np.pi*tail.aspect / ( 2 + np.sqrt(4 + (tail.aspect*beta/etha)**2 * (1 + (np.tan(np.deg2rad(tail.sweep)) / beta) **2 )    )   )
    return CL*np.pi/180.
    
def DATCOM_main(mach,main,width,etha=0.95):
    beta = np.sqrt(1-mach**2)
    CL = 2*np.pi*main.aspect / ( 2 + np.sqrt(4 + (main.aspect*beta/etha)**2 * (1 + (np.tan(np.deg2rad(main.sweep)) / beta) **2 )    )   )
    return (CL*(1+2.15*width/(main.surface/main.chord))* (main.surface-width*main.chord)/main.surface +np.pi/2*width**2/main.surface )*np.pi/180.
    
def downwash(mtv,main,tail):
    r = tail.dist_np/(main.span/2.)
    Kesw = (0.1124 + 0.1265*np.deg2rad(main.sweep)+0.1766*np.deg2rad(main.sweep)**2 )/ r**2 + 0.1024/r +2
    Ke0 = (0.1124)/ r**2 + 0.1024/r +2
    wash = Kesw/Ke0*( r/(r**2 + mtv**2)* 0.4876/(np.sqrt(r**2+0.6319+mtv**2)) \
        +( 1+(r**2/(r**2 + 0.7915+5.0734*mtv**2))**0.3113  )*(1-np.sqrt( (mtv**2/(1+mtv**2)\
        *main.clalpha/(np.pi*main.aspect)) )) )
    return wash*np.pi/180.
    
def calc_xcg(main,tail,canard):
    Mw = 150.
    Sw = main.surface
    xw = main.root/3.    
    MSR = Mw / Sw
    
    Sc = canard.surface
    Mc = MSR*Sc
    xc = canard.dist_np+1./4*main.root
    
    St = tail.surface
    Mt = MSR*St
    xt = tail.dist_np+1./4*main.root
    
    Mland = 100.
    xland = 0.61*main.root
    
    Mprop = 100.
    xprop = -0.01*main.root
    
    Mpl = 10
    xpl = 0.20*main.root
    
    Mbat = 200.
    xbat = 0.3*main.root
    
    Mcomp = 10
    xcomp = 0.17*main.root
    
    Mtotal = Mw+Mc+Mt+Mland+Mpl+Mcomp+Mprop
    
    aft  = 1/Mtotal*(Mw*xw + Mc*xc + Mt+xt + Mland*xland+ Mprop*xprop + Mpl*xpl + Mbat*xbat + Mcomp*xcomp)
    forward = 1/(Mtotal-Mland)*(Mw*xw + Mc*xc + Mt+xt + Mprop*xprop + Mpl*xpl + Mbat*xbat + Mcomp*xcomp)
    return (forward,aft)
    
def est_mass(main,tail,canard):
    Mw = 150.
    Sw = main.surface
    MSR = Mw / Sw
    
    Sc = canard.surface
    Mc = MSR*Sc
    
    St = tail.surface
    Mt = MSR*St
    
    return Mt,Mc
    

def optimize_ratio(xac,canard,main,tail,configuration="t",ran=(0,5),step=0.1):
    #import sympy.plotting
    ratio_range=np.arange(ran[0],ran[1],step)
    mass = []
    for ratio in ratio_range:
        for i in range(0,10):
            stab=sizeStab(xac,canard,main,tail,stabMargin,configuration,ratio)[0]
            cont=sizeControl(xac,canard,main,tail,configuration,ratio)[0]
            
            f1 = sympy.lambdify(stab.get(stab.keys()[0]).free_symbols.pop(),stab.get(stab.keys()[0]))
            f2 = sympy.lambdify(cont.get(cont.keys()[0]).free_symbols.pop(),cont.get(cont.keys()[0]))
        
        
            xcg_left,xcg_right = calc_xcg(main,tail,canard)
            S_l1 = f1(xcg_left)
            S_l2 = f2(xcg_left)
            S_r1 = f1(xcg_right)
            S_r2 = f2(xcg_right)
            #print S_l,S_r
            S = max([S_l1,S_l2,S_r1,S_r2])
            
            if configuration=="t":
                canard.surface= S*ratio
                tail.surface=S
            elif configuration=="c":
                canard.surface= S
                tail.surface=S*ratio
            xcg_left,xcg_right = calc_xcg(main,tail,canard)
        #print est_mass(main,tail,canard)
        mass.append(sum(est_mass(main,tail,canard)))
    #print len(mass)
    #print len(ratio_range)
    plt.plot(ratio_range,mass)
    plt.title("Estimated mass vs Canard-Tail Ratio")
    plt.xlabel(r"$R_{H}$",fontsize=14)
    plt.ylabel(r"Expected Mass [kg]",fontsize=14)
    plt.grid(True)
    ratio_min=ratio_range[np.array(mass).argmin()]
    print ratio_min
    return ratio_min
    
    
def norm_moi(Ixx,Iyy,Izz,Ixz,c,b,mass):
    KX2 = Ixx/(b**2*mass)
    KY2 = Iyy/(c**2*mass)
    KZ2 = Izz/(b**2*mass)
    JXZ = Ixz/(b**2*mass)
    return KX2,KY2,KZ2,JXZ
    
    
if __name__=="__main__":
    
    stabMargin=0.1
    configuration="t"
    ratio=0
    velocity = 100
    mach = velocity/util.scale_a()
    canard,main,tail,vert=dummyWings()
    xac = main.root/4.
    
    xnp = xac + stabMargin*main.root
    
    r_fus=0.6
    t_fus=4./1000
    width_fus=r_fus*2
    l_fus = 5.
    
    mtv = 0.1 # perpendicular distance between zero lift line and tail
    
    rati=np.arange(0,3,0.1)
    mass_r=[]
    
    #for rat in rati:
    for i in range(0,10):
        #tail.clalpha=DATCOM_tail(mach,tail)
        #main.clalpha=DATCOM_main(mach,main,width_fus)
        tail.wash=-downwash(mtv,main,tail)
        xcgs=calc_xcg(main,tail,canard)
        minxcg,min_S,canard.surface,tail.surface = return_sizing(xac,canard,main,tail,configuration,ratio,xcg=xcgs,plot=False)
        canard.chord=(canard.surface/canard.aspect)**0.5
        canard.span=canard.aspect*canard.chord
        tail.chord=(tail.surface/tail.aspect)**0.5
        tail.span=tail.aspect*tail.chord
    #mass_r.append(sum(est_mass(main,tail,canard)))
    #plt.plot(rati,mass_r)
    
    #return_sizing(xac,canard,main,tail,configuration,ratio,xcg=xcgs,plot=True)
    #optimize_ratio(xac,canard,main,tail,ran=[0,5],step=0.001)
    tail.surface=14.
    print("\n")
    #print("DATCOM tail CLalpha: ", DATCOM_tail(mach,tail))
    #print("DATCOM main CLalpha: ", DATCOM_main(mach,main,width_fus))
    print("Emperical Downwash tail: ",-downwash(mtv,main,tail))
    print("xcg with: ",xcgs[0])
    print("xcg without: ",xcgs[1])
    print("Stability: ",sizeStab(xac,canard,main,tail,stabMargin,configuration,ratio))
    print("Control: ",sizeControl(xac,canard,main,tail,configuration,ratio))
    #print("Min S xcg: ",minxcg)
    #print("Minimum S canard: ", canard.surface)
    print("Minimum S tail: ", tail.surface)
    print("Expected Mass: ",sum(est_mass(main,tail,canard)))
    print("dcm/dalpha without: ",cmalpha(canard,main,tail,xac,xcgs[0]))
    print("dcm/dalpha with: ",cmalpha(canard,main,tail,xac,xcgs[1]))
    print("")
    print("Tail Surface: ",tail.surface)
    print("Tail Volume: ",tail.surface*tail.dist_np  /(main.surface*main.chord))
    print("Tail Lift slope: ",tail.clalpha)
    print("Tail wash slope: ",tail.wash)
    print("")
    print("Vert Surface: ",vert.surface)
    print("Vert Volume: ",vert.surface*vert.dist_np  /(main.surface*main.chord))
    print("Vert Lift slope: ",vert.clalpha)
    print("Vert wash slope: ",vert.wash)
    #print("Init. pitch moment derivitive: ",init_v_gust(canard,main,tail,xac,xcg,gust_v,mass,KY2,velocity,density))


    
    Ixx=3500
    Iyy=6000
    Izz=9000
    Ixz=200
    
    mass = 850.
    
    dyn = True
    if dyn:
        co = AircraftStabilityCoeff.coeff()
        
        CL=main.cl
        CD=main.cd
        V0=140.
        rho0=0.45#1.5940
        alpha0=1*np.pi/180.
        
        S=main.surface
        c=main.chord
        b=main.span
        e=main.oswald
        A=main.aspect
        
        
        propInc=0*np.pi/180
        propArm=1.5
        co._steady_conditions(V0,rho0,CD,CL,alpha0)
        co._aircraft_properties(b,c,A,S,e,mass,Ixx,Iyy,Izz,Ixz,xcgs[0],xac,propInc,propArm)
        co._tail(tail)
        co._mainWing(main)
        co._canard(canard)
        co._vert(vert)
        co.deriv()
        
        CZu=0
        Cma=0
        Cmu=0
        co.delta_long(Cma=Cma,CZu=CZu,Cmu=Cmu)
        
        ssS=co.stateSpace(symmetric=True)
        T=np.arange(0,40,0.1)
        u=np.zeros((len(T),2))
        
        upgust = 0#14.#15. #m/s
        alpha = np.arctan(upgust/V0)
        #u[1][0]=alpha
        
        frontgust = 14. #m/s
        du = np.arctan(frontgust/V0)
        #u[0][0]=du
        # uhat, alpha, theta, q*c/V]
        init=[du,alpha,0,0]
        yout,T,xout=control.lsim(ssS,u,T,init)
        ACStab.plot_symmetric(co,yout,T)
        
        
        ssS=co.stateSpace(symmetric=False)
        T=np.arange(0,40,0.01)
        u=np.zeros((len(T),2))
        
        sidegust =14
        
        beta=np.arctan(sidegust/V0)#10*np.pi/180.
        varphi=0
        p=0
        r=0#5*np.pi/180
        init=[beta,varphi,p,r]
        yout,T,xout=control.lsim(ssS,u,T,init)
        #ACStab.plot_asymmetric(co,yout,T)
        
        damp=co.damping(True)
        print "Damping 1: ",damp[0]
        print "Damping 2: ",damp[1]
        
        damp=co.damping(False)
        print "Damping 1: ",damp[0]
        print "Damping 2: ",damp[1]
        print "Damping 3: ",damp[2]
    
        