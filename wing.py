__author__ = 'Stefan'
ji
import numpy as np

b=50.

Cl=10.
Cd=1.
rho=3.
V=40.       #velocity in m/s
q=1./2*rho*V**2
c=[1.,1.,8.,1.,1.,1.,1.] #list of chord lengths with dx steps, from root to tip


def loadcase(Cl,Cd,q,c):
    carray=np.array(c)
    cavg=np.average(carray)
    CL=Cl
    CD=Cd
    n=2*len(c)-1.       #points over all wing
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

    for i in range(len(c)-1):
        x=dx*i
        l=Cl*q*carray[:(i)]       #list of lift/m left to point
        d=Cd*q*carray[:(i)]       #list of drag/m left to point
        V1=L/2.-sum(l)*dx#np.trapz(l,x=None,dx=dx)
        V2=D/2.-sum(d)*dx#np.trapz(d,x=None,dx=dx)
        lnot=Cl*q*carray[(i+1):]    #list of lift/m right to point
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
    return Vy, Vz, Mz, My
    

print loadcase(Cl,Cd,q,c)

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
#
