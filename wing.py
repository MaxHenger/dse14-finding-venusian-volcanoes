__author__ = 'Stefan'

import numpy as np

#b=50.
#
#Cl=0.175
#Cd=0.0075
#rho=3.
#V=40.       #velocity in m/s
#q=1./2*rho*V**2
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
#    for i in range(len(c)-1):
#        W=W+g*rhowing*carray[i]*F*tw*dx  #weight of 1 wing
    for i in range(len(c)):
        l=Cl*q*carray[:i]       #list of lift/m left to point
        #w=g*rhowing*carray[:i]*F*tw   #tw=thickness, F=length of bars in crossection divided by c (constant), weight/m
        w = 0.        
        d=Cd*q*carray[:i]       #list of drag/m left to point
        V1=L/2.-W-sum(l-w)*dx#np.trapz(l,x=None,dx=dx)
        V2=D/2.-sum(d)*dx#np.trapz(d,x=None,dx=dx)
        cloop=carray[i+1:]
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
    """y in lift, z in drag, Mz moment around z, My moment around y"""
    return Vy, Vz, Mz, My, c
  
#A=1.77
#taper=0.55
#F=0.
#tw=0.005
#g=8.
#S=35.
#rhowing=3000.
#c = loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[-1]
#Mz = loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[2]
#Vy = loadcase(Cl,Cd,q,A,taper,S,rhowing,F,tw,g)[0]
