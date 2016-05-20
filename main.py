__author__ = 'Stefan'
import numpy as np
from mayavi import mlab
from itertools import chain
############main##################################333
rhofuselage=4430.#http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MTP642
rhowing=4430.#
sufuselage=1100000000.
suwing=500000.#
V=40.#
rhoair=2.#
Volume=20.#
L=5.#
Loc1=3.1#
Loc2=4.#
taper=0.55
cl=2.#depends on airfoil
cd=0.4#depends on airfoil
A=4.77
S=35.
F=0.004#depends on airfoil,
g=8.87#
SF=1.5
p=0.#   #pressure
E=114*10**9
b=np.sqrt(S*A)
m=18*300*10**(-6)*b/2*rhowing*2   #2 for extra material like bateries, ailerons
Iyyavg=0.00089135 #drag
Ixxavg=2.313445*10**(-5) #lift
ky=78.37*E*Iyyavg/((b/2)**3)
kx=78.37*E*Ixxavg/((b/2)**3)
wny=np.sqrt(ky/m)
wnx=np.sqrt(kx/m)
fny=wny/2/np.pi
fnx=wnx/2/np.pi
print fny, fnx


import fuselage as fs


def structuremass(rhofuselage,rhowing,sufuselage,suwing,V,rhoair,Volume, L, Loc1, Loc2,taper,cl,cd,A,S,F,g,SF,p):#A=aspect ratio
    tf=0.0001    #fuselage thickness
    tw=0.0001    #wingbox thickness
    q=1/2*rhoair*V**2
    R=0.6#fs.fuselagesize(Volume,L)[0]
    #print R
    #rhofuselage=    #density fuselage material
    #rhowing=        #density wingbox material
    deviation=1.
    Lift1=6160.#lift1(M,L,)
    Lift2=0.#M-Lift1#lift2(M,Lift1)
    while deviation>0.0002:


        fuselageloadcase=fs.fuselageloadcase(rhofuselage,g,R,L,tf,Lift1,Lift2,Loc1,Loc2)
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
        deviationfuselage=abs(tf2-tf)
        #deviationwing=np.absolute(tw2-tw)
        deviation=deviationfuselage#+deviationwing
        tf=tf2
        #tw=tw2

    Mf=fs.fuselagemass(tf,L,R,rhofuselage)
    #Mw=wingmass(,,rhowing)#multiply with constant for stuff other than wingbox
    #Mstructure=Mf+Mw
    fuselage=fs.fuselagestress(R,sufuselage,Vyf,Vzf,Mzf,Myf,tf,SF,p,L)
    x=np.array(fuselage[1])
    y=np.array(fuselage[2])
    z=np.array(fuselage[3])
    vMs=np.array(fuselage[4])

    x2=list(chain.from_iterable(x))
    y2=list(chain.from_iterable(y))
    z2=list(chain.from_iterable(z))

    print fs.plot_mayavi([x2,y2,z2,vMs])

    return tf, Mf

    #return Mstructure
print structuremass(rhofuselage,rhowing,sufuselage,suwing,V,rhoair,Volume, L, Loc1, Loc2,taper,cl,cd,A,S,F,g,SF,p)
#print structuremass(1440.,0.,1500000000,.0,3.,3.,20., 5., 1., 4.,0.5,2.,0.5,5.,30.,0.,8.,2.,3000000.)
