__author__ = 'Stefan'
import numpy as np


rho=1000.   #material density
g=10. #grav acceleration
t=0.001
Volume=10.
L=4.

Lift1=50000.        #lift depends on everything
Lift2=10000.
Loc1=0.5
Loc2=3.
print fuselageloadcase(rho,g,R,L,t,Lift1,Lift2,Loc1,Loc2)
#assumptions: fuselage=cilinder+spherical ends, 150 kg of other material in fuselage, only Lift in wing, no weight


############main##################################333
rhofuselage=4000.
rhowing=5000.
Efuselage=500000.
Ewing=500000.
V=40.
Rhoair=2.
Volume=4.
L=5.
Loc1=3.
Loc2=4.
taper=0.6
cl=1.2#depends on airfoil
cd=2.4#depends on airfoil
A=8.
S=20.
F=0.004#depends on airfoil,
g=8.
SF=2.
p=3000000   #pressure






def structuremass(rhofuselage,rhowing,Efuselage,Ewing,V,rhoair,Volume, L, Loc1, Loc2,taper,cl,cd,A,S,F,g,SF):#A=aspect ratio
    tf=0.0001    #fuselage thickness
    tw=0.0001    #wingbox thickness
    q=1/2*rhoair*V**2
    R=fuselagesize(Volume,L)[0]
    #rhofuselage=    #density fuselage material
    #rhowing=        #density wingbox material
    while deviation>0.0002:

        Lift1=lift1(M,L,)
        Lift2=M-Lift1#lift2(M,Lift1)
        fuselageloadcase=fuselageloadcase(rhofuselage,g,R,L,t,Lift1,Lift2,Loc1,Loc2,p)
        wingloadcase=wingloadcase(cl,cd,q,A,taper,S,F,tw,g)#c should depend on amount of lift and aspect ratio, taper ratio, area
        structure=structure()#moments of inertia, E modules etc.,
        Vy=wingloadcase[0]
        Vz=wingloadcase[1]
        Mz=wingloadcase[2]
        My=wingloadcase[3]
        c=wingloadcase[4]
        Vyf=fuselageloadcase[0]
        Vzf=fuselageloadcase[1]
        Mzf=fuselageloadcase[2]
        Myf=fuselageloadcase[3]
        fuselagestresses=fuselagestress(R,sufuselage,Vyf,Vzf,Mzf,Myf,tf,SF,p,L)
        wingstresses=wingstress(I,Ewin,Vy,Vz,Mz,My,structure,tw,c,SF)
        tf2=fuselagestresses
        tw2=wingstresses[0]
        deviationfuselage=np.absolute(tf2-tf)
        deviationwing=np.absolute(tw2-tw)
        deviation=deviationwing+deviationfuselage
        tf=tf2
        tw=tw2

    Mf=fuselagemass(t,L,R,rhofuselage)
    Mw=wingmass(,,rhowing)#multiply with constant for stuff other than wingbox
    Mstructure=Mf+Mw
    fuselage=fuselagestress(R,sufuselage,Vyf,Vzf,Mzf,Myf,tf,SF,p,L)
    x=np.array(fuselage[1])
    y=np.array(fuselage[2])
    z=np.array(fuselage[3])
    vMs=np.array(fuselage[4])
    #dz=fuselage[5]
    
    
    #map(Y.extend,y)
    #map(Z.extend,z)
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    
    #ax.plot_surface(Z,X,Y,rstride=1, cstride=1)#color=vMs
    #plt.legend()
    #plt.show()
    from mayavi import mlab
    from itertools import chain
    X=list(chain.from_iterable(x))
    Y=list(chain.from_iterable(y))
    Z=list(chain.from_iterable(z))
    #print np.shape(Y)
    #print np.shape(Z)
    #print np.shape(vMs)
    def plot_mayavi(structured_data):
        sigma = structured_data[3]
        fig = mlab.figure(bgcolor=(1,1,1),fgcolor=(0,0,0))
        #obj=mlab.contour3d(x, y, z,sigma)
        obj = mlab.points3d(structured_data[0], structured_data[1], structured_data[2],sigma,mode='cube', scale_mode="none",scale_factor=0.3)#mode="cube"
        bar = mlab.colorbar(object=obj, title="Stress (in Pa)", orientation="vertical")
        mlab.view(50,50)
        mlab.roll(5)
        #mlab.savefig(os.path.dirname(os.path.realpath(__file__))+'\\sim_outputs\\plot_mayavi.jpg', size=[5000,3000])
        mlab.show()
    print plot_mayavi([X,Y,Z,vMs])
    return Mstructure
