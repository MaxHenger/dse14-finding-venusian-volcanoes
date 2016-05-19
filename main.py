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
#print fuselageloadcase(rho,g,R,L,t,Lift1,Lift2,Loc1,Loc2)
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




import fuselage as fs


def structuremass(rhofuselage,rhowing,sufuselage,suwing,V,rhoair,Volume, L, Loc1, Loc2,taper,cl,cd,A,S,F,g,SF, p):#A=aspect ratio
    tf=0.0001    #fuselage thickness
    tw=0.0001    #wingbox thickness
    q=1/2*rhoair*V**2
    R=fs.fuselagesize(Volume,L)[0]
    #rhofuselage=    #density fuselage material
    #rhowing=        #density wingbox material
    deviation=1.
    while deviation>0.0002:

        Lift1=30000.#lift1(M,L,)
        Lift2=1000.#M-Lift1#lift2(M,Lift1)
        fuselageloadcase=fs.fuselageloadcase(rhofuselage,g,R,L,t,Lift1,Lift2,Loc1,Loc2,p)
        #wingloadcase=wingloadcase(cl,cd,q,A,taper,S,F,tw,g)#c should depend on amount of lift and aspect ratio, taper ratio, area
        #structure=structure()#moments of inertia, E modules etc.,
        #Vy=wingloadcase[0]
        #Vz=wingloadcase[1]
        #Mz=wingloadcase[2]
        #My=wingloadcase[3]
        #c=wingloadcase[4]
        Vyf=fuselageloadcase[0]
        Vzf=fuselageloadcase[1]
        Mzf=fuselageloadcase[2]
        Myf=fuselageloadcase[3]

        fuselagestresses=fs.fuselagestress(R,sufuselage,Vyf,Vzf,Mzf,Myf,tf,SF,p,L)
        #wingstresses=wingstress(I,Ewin,Vy,Vz,Mz,My,structure,tw,c,SF)
        tf2=fuselagestresses[0]
        #tw2=wingstresses[0]
        deviationfuselage=np.absolute(tf2-tf)
        #deviationwing=np.absolute(tw2-tw)
        deviation=deviationfuselage#+deviationwing
        tf=tf2
        #tw=tw2

    #Mf=fuselagemass(t,L,R,rhofuselage)
    #Mw=wingmass(,,rhowing)#multiply with constant for stuff other than wingbox
    #Mstructure=Mf+Mw
    fuselage=fs.fuselagestress(R,sufuselage,Vyf,Vzf,Mzf,Myf,tf,SF,p,L)
    x=np.array(fuselage[1])
    y=np.array(fuselage[2])
    z=np.array(fuselage[3])
    vMs=np.array(fuselage[4])

    from itertools import chain
    X=list(chain.from_iterable(x))
    Y=list(chain.from_iterable(y))
    Z=list(chain.from_iterable(z))
    #print np.shape(Y)
    #print np.shape(Z)
    #print np.shape(vMs)
    return tf
    print plot_mayavi([X,Y,Z,vMs])
    #return Mstructure

print structuremass(1440.,0.,1500000000,.0,3.,3.,20., 5., 1., 4.,0.5,2,0.5,5.,30,0.,8.,2,8*10**5)
