__author__ = 'Stefan'
import numpy as np


def fuselagesize(V,L):#volume, cilinder length
    L=5.#length cilinder of fuselage
    V=6. #Volume in m3
    r=np.roots([4/3*np.pi,0,2*np.pi*L,-V])
    R=np.real(r[2])#fuselage radius for given length and volume
    tL=L+2*R #total length
    return R, tL


#specificstrength=250000 #Nm/kg
#p=500000  #pressure difference in Pa
#stresses multiplied by t
#Cilinder
#Loopstress=p*R/t
#Longitudinalstress=p*R/2
#sphere
#t=0.001 #initial thickness
#Spherestress=p*R/2 #loop or longitudinal
#VMcilinder=np.sqrt(1./2*(Longitudinalstress**2+Loopstress**2+(Longitudinalstress-Loopstress)**2))
#VMsphere=Spherestress #von mises stress of sphere is equal to spherestress


rho=1000.   #material density
g=10. #grav acceleration

def fuselageloadcase(rho,g,R,L,t,Lift1,Lift2,Loc1,Loc2):#lift1=lift main wing-weight main wing, Loc1=location of main wing, Lift2=lift second  wing
    tL=2*R+L
    #z pointed towards the back
    #y positive down
    #discretization
    n=10
    dz=tL/n
    #weight distribution+lift
    wd=[]   #weight in N, positive downwards
    for i in range(n):
        z=i*dz
        if z<R:
            wdl=0.1
        elif z<(R+L):
            wdl=t*2*np.pi*R*g*rho
        else:# z<(2*R+L):
            wdl=0.1
        if z<Loc1<(z+dz):
            wdl=wdl-Lift1
        if z<Loc2<(z+dz):
            wdl=wdl-Lift2
        wd.append(wdl)
    Vx=[]
    Vy=[]
    My=[]
    Mx=[]
    V1=0.
    for i in range(n):
        z=i*dz
        V1=V1+dz*wd[i]
        V2=0.
        M1=0.
        M2=0.
        for i in range(n):
            M1=M1+z*dz*wd[i]
            M2=M2+0.

        Vx.append(V1)
        Vy.append(V2)
        My.append(M1)
        Mx.append(M2)
    return Vx,Vy, My, Mx
rho=1000.   #material density
g=10. #grav acceleration
t=0.001
Volume=10.
L=4.
R=fuselagesize(Volume,L)[0]
Lift1=50000.        #lift depends on everything
Lift2=10000.
Loc1=0.5
Loc2=3.
print fuselageloadcase(rho,g,R,L,t,Lift1,Lift2,Loc1,Loc2)


t=0.001
rho=3000.
cl=
cd=
q=
rho=
while deviation<0.001:
    M=mass(t,rho,)
    Lift1=lift1(M,L,)
    Lift2=M-Lift1#lift2(M,Lift1)
    fuselageloadcase=fuselageloadcase(rho,g,R,L,t,Lift1,Lift2,Loc1,Loc2)
    wingloadcase=wingloadcase()
    structure=structure(cl,cd,q,A)#moments of inertia, E modules etc., c should depend on amount of lift and aspect ratio or depend on area
    stresses=stress(fuselageloadcase,wingloadcase,structure)
    t=
