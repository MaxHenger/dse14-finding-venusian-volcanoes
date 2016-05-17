__author__ = 'Stefan'
import numpy as np



specificstrength=250000 #Nm/kg
p=500000  #pressure difference in Pa
t=0.001 #initial thickness
L=5.#length cilinder of fuselage

V=6. #Volume in m3
r=np.roots([4/3*np.pi,0,2*np.pi*L,-V])
print r
R=np.real(r[2])#fuselage radius for given length and volume
tL=L+2*R #total length
#stresses multiplied by t
#Cilinder
Loopstress=p*R/t
Longitudinalstress=p*R/2
#sphere
print Longitudinalstress
print Loopstress
Spherestress=p*R/2 #loop or longitudinal
VMcilinder=np.sqrt(1./2*(Longitudinalstress**2+Loopstress**2+(Longitudinalstress-Loopstress)**2))
VMsphere=Spherestress #von mises stress of sphere is equal to spherestress
print VMcilinder
print VMsphere