# -*- coding: utf-8 -*-
"""
Created on Mon May 30 10:40:22 2016

Calculates the coefficients needed to calculate / simulate the aircraft.
Calculates the matrices

@author: Julius
"""

import numpy as np
import control

class coeff:
    def __init__(self):
        self.ACpropCheck=False
        self.SSFlightCheck=False
    
    def _steady_conditions(self,V0,rho0,CD,CL,alpha0,grav=8.75):
        self.V0=V0
        self.rho0=rho0
        self.CL=CL
        self.CD=CD
        self.Tc=CD # thrust = drag
        self.alpha0=alpha0
        self.gamma0=0
        self.g=grav
        self.SSFlightCheck=True
    
    def _aircraft_properties(self,span,chord,aspect,SurfaceArea,oswald,mass,Ixx,Iyy,Izz,Ixz,xcg,xac,propInc,propArm):
        self.b=span
        self.c=chord
        self.A=aspect
        self.S=SurfaceArea
        self.mass=mass
        self.e=oswald
        self.xcg=xcg
        self.xac=xac
        self.Ixx=Ixx
        self.Iyy=Iyy
        self.Izz=Izz
        self.Ixz=Ixz
        self.propInc=propInc
        self.propArm=propArm
        self.ACpropCheck=True      
        
    def _mainWing(self,wing):
        self.main=wing
    
    def _canard(self,wing):
        self.canard=wing
    
    def _tail(self,wing):
        self.tail=wing
    
    def _vert(self,wing):
        self.vert=wing
    
    def MatrixA(self,symmetric=True):
        if symmetric:
            C1, C2, C3 =  self._symmetricMatrices()
        else:
            C1, C2, C3 = self._asymmetricMatrices()
        self.MA = C1**-1 * C2
        return self.MA
        
    def MatrixB(self,symmetric=True):
        if symmetric:
            C1, C2, C3 =  self._symmetricMatrices()
        else:
            C1, C2, C3 = self._asymmetricMatrices()
        self.MB = C1**-1 * C3
        return self.MB
        
    def _symmetricMatrices(self):
        SC1 = np.matrix([[-2.*self.muc*self.c/self.V0, 0,                     0,       0],
                         [0,             (self.CZadot-2.*self.muc)*self.c/self.V0,   0,       0],
                         [0,              0,                     -self.c/self.V0,    0],
                         [0,             self.Cmadot*self.c/self.V0,           0,       -2.*self.muc*self.KY2*(self.c/self.V0)]])
        
        SC2 = np.matrix([[ -self.CXu , -self.CXa , -self.CZ0, -self.CXq  ],
                         [ -self.CZu , -self.CZa , self.CX0 , -(self.CZq+2*self.muc)],
                         [ 0 ,0 ,0, -1],
                         [ -self.Cmu , -self.Cma, 0, -self.Cmq]])
                                
        SC3 = np.matrix([[-self.CXde,-self.CXde],
                         [-self.CZde,-self.CZde],
                         [0.,0],
                         [-self.Cmde,-self.Cmde]])
        return SC1, SC2, SC3
    
    def _asymmetricMatrices(self):
        AC1=np.matrix([[(self.CYbdot-2.*self.mub)*self.b/self.V0,      0,                0,              0],
                  [0,                            -0.5*self.b/self.V0,         0,              0],
                  [0,                            0,               -4.*self.mub*self.KX2*(self.b/self.V0),    4.*self.mub*self.KXZ*(self.b/self.V0)],
                  [self.Cnbdot*self.b/self.V0,                 0,                4.*self.mub*self.KXZ*(self.b/self.V0),    -4.*self.mub*self.KZ2*(self.b/self.V0)]])
        
        AC2=np.matrix([[ -self.CYb    , -self.CL      , -self.CYp     , -(self.CYr - 4*self.mub) ],
                   [ 0          , 0       ,  -1           , 0                    ],
                   [ -self.Clb        , 0       , -self.Clp       , -self.Clr            ],
                   [ -self.Cnb        , 0       , -self.Cnp       , -self.Cnr           ]])
            
        AC3=np.matrix([[ -self.CYda   , -self.CYdr    ],
                   [ 0      , 0       ],
                   [ -self.Clda   , -self.Cldr    ],
                   [ -self.Cnda   , -self.Cndr    ]])
    
        return AC1, AC2, AC3
    
    def stateSpace(self,symmetric=True):
    
        MA = self.MatrixA(symmetric)
        MB = self.MatrixB(symmetric)
        MC = np.eye(4)
        MD = np.zeros(MB.shape)
        
        return control.ss(MA,MB,MC,MD)    
    
    def damping(self,symmetric=True):
        ssS=self.stateSpace(symmetric)
        poles=ssS.pole()
        print poles
        if symmetric:
            return -poles[0].real/abs(poles[0]), -poles[2].real/abs(poles[2])
        else:
            return -poles[0].real/abs(poles[0]), -poles[1].real/abs(poles[1]), -poles[3].real/abs(poles[3])
    
    def _mass_coeff(self):
        self.W=self.mass*self.g
        self.muc  = self.mass / (self.rho0 * self.S * self.c)
        self.mub  = self.mass / (self.rho0 * self.S * self.b)
        self.KY2 = self.Iyy/(self.mass*self.c**2)
        self.KX2 = self.Ixx/(self.mass*self.b**2)
        self.KZ2 = self.Izz/(self.mass*self.b**2)
        self.KXZ = self.Ixz/(self.mass*self.b**2)
        
    def dTdV(self): # not used
        return -3*self.Tc/self.V0
        
    def dCDdV(self): # not used
        return (-4/self.V0**5*self.W**2/((0.5*self.rho0*self.S)**2 *np.pi*self.A*self.e) )
        
    def dCLdV(self): # not used
        return -2/self.V0**3 * self.W/(0.5*self.rho0*self.S)
        
    def dCDdT(self):
        return 0.00
    
    def dCmdT(self):
        return -self.propInc*self.propArm
        
    def deriv(self):
        self._mass_coeff()
        self._deriv_0()
        self._deriv_u()
        self._deriv_a()
        self._deriv_adot()
        self._deriv_q()
        self._deriv_de()
        self._derive_b()
        self._deriv_bdot()
        self._deriv_p()
        self._deriv_r()
        self._deriv_da()
        self._deriv_dr()
        
    def delta_long(self,CX0=0,CZ0=0,Cm0=0,CXu=0,CZu=0,Cmu=0,CXa=0,CZa=0,Cma=0,CXadot=0,CZadot=0,Cmadot=0,CXq=0,CZq=0,Cmq=0,CXde=0,CZde=0,Cmde=0):
        self.CX0+=CX0
        self.CZ0+=CZ0
        self.Cm0+=Cm0
        self.CXu+=CXu
        self.CZu+=CZu
        self.Cmu+=Cmu
        self.CXa+=CXa
        self.CZa+=CZa
        self.Cma+=Cma
        self.CXadot+=CXadot
        self.CZadot+=CZadot
        self.Cmadot+=Cmadot
        self.CXq+=CXq
        self.CZq+=CZq
        self.Cmq+=Cmq
        self.CXde+=CXde
        self.CZde+=CZde
        self.Cmde+=Cmde
        
    def _deriv_0(self):
        self.CX0 = self.W/(0.5*self.rho0*self.V0**2*self.S) * np.sin(self.gamma0)
        self.CZ0 = -self.W/(0.5*self.rho0*self.V0**2*self.S) * np.cos(self.gamma0)
        self.Cm0 = 0
        print "CX0: ",self.CX0
        print "CZ0: ",self.CZ0
        print "Cm0: ",self.Cm0
    def _deriv_u(self):
        self.CXu = 2*self.CL*np.tan(self.gamma0)-3*self.CD*(1-self.dCDdT())
        self.CZu = -2*self.CL+self.CD*(-(self.alpha0+self.propInc) + 3*self.dCDdT() )
        self.Cmu = -3*self.CD*self.dCmdT()
        print "CXu: ",self.CXu
        print "CZu: ",self.CZu
        print "Cmu: ",self.Cmu
    def _deriv_a(self):
        # CXa is negative due to the sign of the CL
        self.CXa = -self.CL*(1-2*self.main.clalpha*(180/np.pi)/(np.pi*self.main.aspect*self.main.oswald))
        self.CZa = -self.main.clalpha\
        -self.tail.clalpha*(1+self.tail.wash)*self.tail.VelFrac**2*self.tail.surface/self.main.surface\
        -self.canard.clalpha*(1-self.canard.wash)*self.canard.VelFrac**2*self.canard.surface/self.main.surface
        self.Cma = self.main.clalpha*(self.xcg-self.xac)/self.main.chord \
         - self.tail.clalpha*self.tail.dist_np*self.tail.surface/(self.main.surface*self.main.chord)*self.tail.VelFrac**2 \
         - self.canard.clalpha*self.canard.dist_np*self.canard.surface/(self.main.surface*self.main.chord)*self.canard.VelFrac**2
        print "CXa: ",self.CXa
        print "CZa: ",self.CZa
        print "Cma: ",self.Cma
         
    def _deriv_q(self):
        self.CXq = 0
        # for CZq,Cmq it is positive in the formulea below due to signs of the cl
        self.CZq = +2*self.tail.cl*self.tail.VelFrac**2*self.tail.surface*self.tail.dist_np/(self.main.surface*self.main.chord)\
                   +2*self.canard.cl*self.canard.VelFrac**2*self.canard.surface*self.canard.dist_np/(self.main.surface*self.main.chord)
        self.Cmq = +1.15*self.tail.cl*self.tail.VelFrac**2*self.tail.surface*self.tail.dist_np**2/(self.main.surface*self.main.chord**2)\
                   +1.15*self.canard.cl*self.canard.VelFrac**2*self.canard.surface*self.canard.dist_np**2/(self.main.surface*self.main.chord**2)
        print "CXq: ",self.CXq
        print "CZq: ",self.CZq
        print "Cmq: ",self.Cmq
    def _deriv_adot(self):
        self.CXadot = 0
        self.CZadot = -self.tail.cl*self.tail.VelFrac**2*self.tail.wash*self.tail.surface*self.tail.dist_np/(self.main.surface*self.main.chord)\
                      -self.canard.cl*self.canard.VelFrac**2*self.canard.wash*self.canard.surface*self.canard.dist_np/(self.main.surface*self.main.chord)
        # sign change due to convention
        self.Cmadot = +self.tail.cl*self.tail.VelFrac**2*self.tail.wash*self.tail.surface*self.tail.dist_np**2/(self.main.surface*self.main.chord**2)\
                      +self.canard.cl*self.canard.VelFrac**2*self.canard.wash*self.canard.surface*self.canard.dist_np**2/(self.main.surface*self.main.chord**2)
        print "CXadot: ",self.CXadot
        print "CZadot: ",self.CZadot
        print "Cmadot: ",self.Cmadot
    def _deriv_de(self):
        self.CXde = 0
        self.CZde = 0#+self.tail.clde*self.tail.VelFrac**2*self.tail.surface/self.main.surface\
                    #+self.canard.clde*self.tail.VelFrac**2*self.canard.surface/self.main.surface
        self.Cmde = -1.08556793814#+self.tail.clde*self.tail.VelFrac**2*self.tail.surface*self.tail.dist_np/(self.main.surface*self.main.chord)\
                    #+self.canard.clde*self.tail.VelFrac**2*self.canard.surface*self.canard.dist_np/(self.main.surface*self.main.chord)
        print "CXde: ",self.CXde
        print "CZde: ",self.CZde
        print "Cmde: ",self.Cmde
        
    # LATERAL DERIVATIVES
    def _derive_b(self):
        self.CYb_v=-self.vert.clalpha*(1+self.vert.wash)*self.vert.VelFrac**2*self.vert.surface/self.main.surface
        self.CYb_prop=0
        self.CYb_fus=0
        self.CYb=self.CYb_v+self.CYb_prop+self.CYb_fus    
        
        self.Clb_w=-1./4* self.main.cl*np.sin(2*np.deg2rad(self.main.sweep))
        self.Clb_v=self.CYb_v*(self.vert.dist_z/self.main.span*np.cos(self.alpha0)-self.vert.dist_np/self.main.span*np.sin(self.alpha0))
        self.Clb=self.Clb_v+self.Clb_w
    
        #a = 0.105*self.main.aspect/(self.main.aspect+2)
        #self.Clb_di=-a/6.*(1+2*self.main.taper)/(1+self.main.taper)*self.main.dihedral
        #self.Clb_sw=-1./3*(1+2*self.main.taper)/(1+self.main.taper)*self.main.cl*np.tan(np.deg2rad(self.main.sweep))
        #self.Clb=self.Clb_di+self.Clb_sw          
        
        self.Cnb_v=self.CYb_v*(self.vert.dist_z/self.main.span*np.sin(self.alpha0)-self.vert.dist_np/self.main.span*np.cos(self.alpha0))
        self.Cnb_prop=0
        self.Cnb_w=0
        self.Cnb=self.Cnb_v+self.Cnb_prop+self.Cnb_w
        
        print "CYb: ",self.CYb
        print "Clb: ",self.Clb
        print "Cnb: ",self.Cnb
        
    def _deriv_bdot(self):
        self.CYbdot=0
        self.Cnbdot=0
        
        print "CYbdot: ",self.CYbdot
        print "Cnbdot: ",self.Cnbdot
    
    def _deriv_p(self):
        self.CYp = -2*self.vert.clalpha*self.vert.dist_z/self.main.span*self.vert.VelFrac**2* self.vert.surface/self.main.surface
        self.Clp_v = -2*self.vert.clalpha*(self.vert.dist_z/self.main.span)**2*self.vert.VelFrac**2*self.vert.surface/self.main.surface
        self.Clp = self.Clp_v -0.32 #-0.35 # emperically found cause fuck this shit, page 220 of reader
        self.Cnp = -2*self.vert.clalpha*self.vert.dist_z*self.vert.dist_np/self.main.span**2*self.vert.VelFrac**2*self.vert.surface/self.main.surface
        print "CYp: ",self.CYp
        print "Clp: ",self.Clp
        print "Cnp: ",self.Cnp
    def _deriv_r(self):
        
        self.CYr_v = 2*self.vert.clalpha*self.vert.VelFrac**2*self.vert.dist_np/self.main.span*self.vert.surface/self.main.surface
        self.CYr_prop=0
        self.CYr_fus=0
        self.CYr = self.CYr_prop+self.CYr_v+self.CYr_fus
        
        self.Clr_v = self.CYr_v*(self.vert.dist_z*np.cos(self.alpha0)-self.vert.dist_np*np.sin(self.alpha0))/self.main.span 
        self.Clr_fus = 0
        self.Clr=self.Clr_v+self.Clr_fus
        
        self.Cnr_v = -self.CYr_v*self.vert.dist_np/self.main.span
        self.Cnr = self.Cnr_v
        
        print "CYr: ",self.CYr
        print "Clr: ",self.Clr
        print "Cnr: ",self.Cnr
    def _deriv_da(self):
        self.CYda=0
        self.Clda=0.19620358663 #-0.23
        self.Cnda=0
        
        print "CYda: ",self.CYda
        print "Clda: ",self.Clda
        print "Cnda: ",self.Cnda
        
    def _deriv_dr(self):
        self.CYdr= self.vert.clde*self.vert.VelFrac**2*self.vert.surface/self.main.surface
        self.Cldr= self.CYdr*(self.vert.dist_z*np.cos(self.alpha0)-self.vert.dist_np*np.sin(self.alpha0))/self.main.span 
        self.Cndr= -self.CYdr*self.vert.dist_np/self.main.span
        
        print "CYdr: ",self.CYdr
        print "Cldr: ",self.Cldr
        print "Cndr: ",self.Cndr
        
if __name__ == "__main__":
    import stability
    canard,main,tail,vert=stability.dummyWings()
    co=coeff()
    V0=100.
    rho0=1.5940
    CL=main.cl
    CD=main.cd
    alpha0=1*np.pi/180.
    S=main.surface
    c=main.chord
    b=main.span
    e=main.oswald
    A=main.aspect
    xcg = 1.367
    xac = 1.23
    m=660
    Ixx=1000
    Iyy=5000
    Izz=1000
    Ixz=0
    propInc=1*np.pi/180
    propArm=4
    co._steady_conditions(V0,rho0,CD,CL,alpha0)
    co._aircraft_properties(b,c,A,S,e,m,Ixx,Iyy,Izz,Ixz,xcg,xac,propInc,propArm)
    co._tail(tail)
    co._mainWing(main)
    co._canard(canard)
    co._vert(vert)
    co.deriv()
    
    damp=co.damping()
    print "Damping 1: ",damp[0]
    print "Damping 2: ",damp[1]
    MA=co.MatrixA()
    MB=co.MatrixB()
    MC = np.eye(4)
    MD = np.zeros(MB.shape)
    
    SS=control.ss(MA,MB,MC,MD)
    
"""
CX0    = W * np.sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu    = -0.02792
CXa    = -0.47966
CXadot = +0.08330
CXq    = -0.28170
CXde   = -0.03728

CZ0    = -W * np.cos(th0) / (0.5 * rho * V0 ** 2 * S)
CZu    = -0.37616
CZa    = -5.74340
CZadot = -0.00350
CZq    = -5.66290
CZde   = -0.69612

Cmu    = +0.06990
Cmadot = +0.17800
Cmq    = -8.79415

CYb    = -0.7500
CYbdot =  0.     
CYp    = -0.0304
CYr    = +0.8495
CYda   = -0.0400
CYdr   = +0.2300

Clb    = -0.10260
Clp    = -0.71085
Clr    = +0.23760
Clda   = -0.23088
Cldr   = +0.03440

Cnb    =  +0.1348
Cnbdot =   0     
Cnp    =  -0.0602
Cnr    =  -0.2061
Cnda   =  -0.0120
Cndr   =  -0.0939
"""