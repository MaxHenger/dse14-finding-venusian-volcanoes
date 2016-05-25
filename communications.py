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
LF_atten = attenuation loss factor <=1 [-]
Gain_rec = receiver gain [-]
freq = frequency [Hz]
distance = distance between transmitter and receiver [m]
pointoff_rec = pointing offset angle receiver [deg]
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
        G_tx = 20*np.log10(D_transm) + 20*np.log10(freq/(1.*10**9)) + 17.8
        halfpowerangle_transm = 21./((freq/(1.*10**9))*D_transm)
        L_ptx = -12 * (pointoff_transm/halfpowerangle_transm)**2 #pointing loss
        
        # Space part        
        L_a = 10*np.log10(LF_atten) #attenuation loss
        c = 299792458.
        wavelength = c/freq  
        L_s = 10*np.log10((wavelength/(4*np.pi*distance))**2)  
        
        # Receiver
        halfpowerangle_rec = 21./((freq/(1.*10**9))*D_rec)
        G_rx = 20*np.log10(D_rec) + 20*np.log10(freq/(1.*10**9)) + 17.8
        L_prx = -12 * (pointoff_rec/halfpowerangle_rec)**2
        L_rx = 10*np.log10(LF_rec)
        
        Boltzmann = 1.38064852*(10**(-23))        
        L_Boltz = 10.*np.log10(1./Boltzmann)
            
            
        SignToNoise = P + L_l + G_tx + L_ptx + L_a + G_rx + L_s + L_prx + L_rx + L_Boltz - 10*np.log10(DataRate) - 10*np.log10(T_syst)
        
        return SignToNoise
        
        
if __name__ == "__main__":
    
 
    
    P_transm = 5.
    LF_transm = 0.8
    D_transm = 0.25
    eff_transm = 0.55
    pointoff_transm = 0.5
    LF_atten = 0.95
    
    freq = 2.20*10**9
    distance = 1.97*10**6
    
    pointoff_rec = 0.
    D_rec = 15.
    eff_rec = 0.55       
    LF_rec = 0.7
    DataRate = 2.27*10**9
    T_syst = 135.

    SNratio = communication(P_transm,LF_transm,D_transm,eff_transm,pointoff_transm,LF_atten,\
                  freq,distance,pointoff_rec,D_rec,eff_rec,LF_rec,DataRate,T_syst)
    print SNratio
    
    plt.plot(D_transm,SNratio)
    plt.show()