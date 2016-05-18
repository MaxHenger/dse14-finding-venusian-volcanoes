__author__ = 'Stefan'

import numpy as np


Cl=4.
Cd=1.
rho=1.
V=40.       #velocity in m/s
q=1./2*rho*V**2
#c=[1.,1.,8.,1.,1.,1.,1.] #list of chord lengths with dx steps, from root to tip


def loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g):
    n=101.#points over all wing
    b=np.sqrt(A*S)
    cavg=S/b
    cr=2*cavg/(1+taper)
    ct=cr*taper
    c=np.arange(cr,ct-(cr-ct)/(n/2+1),-(cr-ct)/((n-1)/2))
    carray=np.array(c)

    CL=Cl
    CD=Cd
    #n=2*len(c)-1.       #points over all wing
    dx=b/(n-1)
    #for i in range(len(c)):
    #    CL=CL+2*Cl*c[i]/dx
    #    CD=CD+2*Cd*c[i]*dx
    L=CL*q*b*cavg       #total lift
    D=CD*q*b*cavg       #total drag
    Vy=[]       #shear due to lift
    Vz=[]       #shear due to drag
    Mz=[]       #moment due to lift
    My=[]       #moment due to drag
    W=0.
    for i in range(len(c)-1):
        W=W+g*rhowing*carray[i]*F*tw*dx  #weight of 1 wing
    for i in range(len(c)-1):
        l=Cl*q*carray[:(i)]       #list of lift/m left to point
        w=g*rhowing*carray[:i]*F*tw   #tw=thickness, F=length of bars in crossection divided by c (constant), weight/m
        d=Cd*q*carray[:(i)]       #list of drag/m left to point
        V1=L/2.-W-sum(l-w)*dx#np.trapz(l,x=None,dx=dx)
        V2=D/2.-sum(d)*dx#np.trapz(d,x=None,dx=dx)
        cloop=carray[(i+1):]
        Mz0=0.
        My0=0.
        for i in range(len(cloop)):
            x=dx*i
            Mz0=Mz0+dx*x*Cl*q*cloop[i]      #moment due to lift
            My0=My0+dx*x*Cd*q*cloop[i]      #moment due to drag
        M1=Mz0
        M2=My0
        Vy.append(V1)
        Vz.append(V2)
        Mz.append(M1)
        My.append(M2)
    print L
    return Vy, Vz, Mz, My, c
    
A=10.
taper=0.5
F=0.006
tw=0.003
g=8.
S=10.
rhowing=3000.
print loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[2]

#def stress(Vy,Vz,Mz,My,T):
#    t1=     #top and bottom
#    t2=     #sides
#
#    h=      #dependend on c
#    w=      #dependend on c
#    Ixx=
#    Iyy=
#    Q=
#    #shear flow
#    q1=     #top and bottom
#    q2=     #sides
#    #shear stress
#    tau1=max(q1)*t1     #top and bottom
#    tau2=max(q2)q*t2        #sides
#    sigmaxz=    #normal stress due to Mz
#    sigmaxy=    #normal stress  du to My
#    sigmax1= #maximum axial stress top and bottom
#    sigmax2=    #maximum axial stress sides
#    #von Mises
#    vMs1=np.sqrt(sigmax1**2+3*tau1)  #top and bottom
#    vMs2=np.sqrt(sigmax2**2+3*tau2)  #sides
