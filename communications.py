"""
Created on Thu May 12 22:06:19 2016

@author: Mathijs

First order communications analysis Venus aircraft assuming parabolic antenna on both ends of the communication


# Inputs
P_transm = transmitter power [W]
LF_transm = transmitter loss factor <=1 [-]
D_transm = transmitter antenna diameter [m]
eff_transm = transmitter antenna efficiency, 0.55 for parabolic [-] 
pointoff_transm = pointing offset angle transmitter [deg]
halfpowerangle_rec = half power angle transmitter [deg]
LF_atten = attenuation loss factor <=1 [-]
Gain_rec = receiver gain [-]
freq = frequency [Hz]
distance = distance between transmitter and receiver [m]
pointoff_rec = pointing offset angle receiver [deg]
halfpowerangle_rec = half power angle receiver [deg]
LF_rec = receiver loss factor [-]
DataRate = data rate [bits/sec]
T_syst = system temperature [K]

Output
SignToNoise = Signal-to-Noise ratio [dB]

"""

import numpy as np
import matplotlib.pyplot as plt

def communication(P_transm,LF_transm,D_transm,eff_transm,pointoff_transm,LF_atten,\
                  freq,distance,pointoff_rec,D_rec,eff_rec,LF_rec,DataRate,T_syst):
                      
        # Transmitter end              
        P = 10*np.log10(P_transm)
        L_l = 10*np.log10(LF_transm)
        G_t = 20*np.log10(D_transm) + 20*np.log10(freq) + 17.8
        halfpowerangle_transm = 21./((freq/(1.*10**9))*D_transm)
        L_pt = -12 * (pointoff_transm/halfpowerangle_transm)**2
        
        # Space part        
        L_a = 10*np.log10(LF_atten)
        c = 299792458.
        wavelength = c/freq  
        L_s = 10*np.log10((wavelength/(4*np.pi*distance))**2)  
        
        # Receiver
        halfpowerangle_rec = 21./((freq/(1.*10**9))*D_rec)
        G_r = 20*np.log10(D_rec) + 20*np.log10(freq) + 17.8
        L_pr = -12 * (pointoff_rec/halfpowerangle_rec)**2
        L_r = 10*np.log10(LF_rec)
        
        Boltzmann = 1.38064852*(10**(-23))        
        L_Boltz = 10*np.log10(Boltzmann)
        print P
        print L_l
        print G_t
        print L_pt
        
        SignToNoise = P + L_l + G_t + L_pt + L_a + G_r + L_s + L_pr + L_r + L_Boltz - 10*np.log10(DataRate) - 10*np.log10(T_syst)
        
        return SignToNoise
        
        
if __name__ == "__main__":
    
 
    
    P_transm = 10.
    LF_transm = 0.9
    D_transm = 0.8 
    eff_transm = 0.55
    pointoff_transm = 4.
    LF_atten = 0.95
    
    freq = 8.40*10**9
    distance = 1.*10**6
    
    pointoff_rec = 1.
    D_rec = 2.5
    eff_rec = 0.55       
    LF_rec = 0.8
    DataRate = 228000
    T_syst = 600.

    D_transm = diam
    SNratio = communication(P_transm,LF_transm,D_transm,eff_transm,pointoff_transm,LF_atten,\
                  freq,distance,pointoff_rec,D_rec,eff_rec,LF_rec,DataRate,T_syst)
    print SNratio
    
    plt.plot(D_transm,SNratio)
    plt.show()