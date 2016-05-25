# -*- coding: utf-8 -*-
"""
Created on Wed May 18 15:34:27 2016

@author: Julius
"""
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

def dummyWings():
    canard = __dummy_wing__()
    #canard.surface=2.
    canard.dist_np=-2
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
    #tail.surface=0
    tail.dist_np=6
    tail.cl = -0.23
    tail.clalpha=0.0705
    tail.wash = -0.1
    tail.VelFrac = 0.8
    return canard,main,tail


def sizeStab(xac,canard,main,tail,stabMargin=0.1,configuration="both",ratio=0):
    """Returns area of either canard, tail, or coefficients of the equation ret= both, c, t
       Input: Canard, Wing, Tail location and aerodynamic properties
       Ouput: Canard, Tail surface Area"""
    
    import sympy
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
    import sympy
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
        return sympy.solve(eq.subs(StA,ratio*ScA))
    elif configuration=="t":
        return sympy.solve(eq.subs(ScA,ratio*StA))
    else:
        raise ValueError("Configuration Unknown")
        
def plot_sizing(xac,canard,main,tail,configuration="t",ratio=0):
    """plots the required surface area as function of x cg
    """
    import sympy.plotting
    stab=sizeStab(xac,canard,main,tail,stabMargin,configuration,ratio)[0]
    cont=sizeControl(xac,canard,main,tail,configuration,ratio)[0]
    
    rang=(stab.get(stab.keys()[0]).free_symbols.pop(),-2,10)
    p1 = sympy.plotting.plot(stab.get(stab.keys()[0]),rang,show=False, line_color='b',ylabel="Surface Area")
    p2 = sympy.plotting.plot(cont.get(cont.keys()[0]),rang,show=False, line_color='r')
    p1.extend(p2)
    p1.show()

if __name__=="__main__":
    
    xac = 1.23
    stabMargin=0.1
    configuration="t"
    ratio=0
    canard,main,tail=dummyWings()
    
    
    print("\n")
    print("Stability: ",sizeStab(xac,canard,main,tail,stabMargin,configuration,ratio))
    print("Control: ",sizeControl(xac,canard,main,tail,configuration,ratio))
    
    plot_sizing(xac,canard,main,tail,configuration,ratio)