"""
Created on Thu May 12 22:06:19 2016

@author: Mathijs

First order communications analysis Venus aircraft assuming parabolic antenna on both ends of the communication


# Inputs
P_transm = transmitter power [W]
LF_transm = transmitter loss factor <=1 [-]
D_transm = transmitter antenna diameter/horn size/helical diameter [m]
    dish = antenna diameter
    horn = horn size
    helical = diameter of helix
pointoff_transm = pointing offset angle transmitter [deg]
LF_atten = attenuation loss factor <=1 [-]
Gain_rec = receiver gain [-]
freq = frequency [Hz]
distance = distance between transmitter and receiver [m]
pointoff_rec = pointing offset angle receiver [deg]
D_rec = Receiver antenna diameter/horn size/helical diameter [m]
    dish = antenna diameter
    horn = horn size
    helical = diameter of helix
LF_rec = receiver loss factor [-]
DataRate = data rate [bits/sec]
T_syst = system temperature [K]
anttype_TX = dish/horn/helical defines function what antenna is used for transmitter
anttype_RX = dish/horn/helical defines function what antenna is used for receiver
L_helical_TX= length of helix if helix antenna is used for transmitter [m]
L_helical_RX= length of helix if helix antenna is used for receiver [m]



Output
SignToNoise = Signal-to-Noise ratio [dB]



note if using a helical antenna for receiving end set the transmitting helical length to 0
"""

import numpy as np
import matplotlib.pyplot as plt

def communication(P_transm,LF_transm,D_transm,pointoff_transm,LF_atten,\
                  freq,distance,pointoff_rec,D_rec,LF_rec,DataRate,T_syst, anttype_TX, anttype_RX,L_helical_TX=0, L_helical_RX=0):
         

        c = 299792458.
        wavelength = c/freq   

        freq = freq/(1.*10**9)          
        # Transmitter end              
        P = 10*np.log10(P_transm)
        L_l = 10*np.log10(LF_transm)
        if anttype_TX == "dish":
            G_tx = 20*np.log10(D_transm) + 20*np.log10(freq) + 17.8
            halfpowerangle_transm = 21./(freq*D_transm)            
        elif anttype_TX == "horn":
            G_tx = 20*np.log10(np.pi*D_transm/wavelength)-2.8
            halfpowerangle_transm = 225./(np.pi*D_transm/wavelength)            
        elif anttype_TX == "helical":
            G_tx = 10*np.log10(np.pi**2*D_transm**2*L_helical_TX/wavelength**3)+10.3
            halfpowerangle_transm = 52./np.sqrt(np.pi**2*D_transm**2*L_helical_TX/wavelength**3)
        else:
            raise ValueError("wrong anttype input for TX")
        L_ptx = -12. * (pointoff_transm/halfpowerangle_transm)**2 #pointing loss
        
        
        
        # Space part        
        L_a = 10*np.log10(LF_atten) #attenuation loss          
        L_s = 10*np.log10((wavelength/(4*np.pi*distance))**2)  
        
        # Receiver        
        if anttype_RX == "dish":
            G_rx = 20*np.log10(D_rec) + 20*np.log10(freq) + 17.8
            halfpowerangle_rec = 21./(freq*D_rec)            
        elif anttype_RX == "horn":
            G_rx = 20*np.log10(np.pi*D_rec/wavelength)-2.8
            halfpowerangle_rec = 225./(np.pi*D_rec/wavelength)            
        elif anttype_RX == "helical":
            G_rx = 10*np.log10(np.pi**2*D_rec**2*L_helical_RX/wavelength**3)+10.3
            halfpowerangle_rec = 52./np.sqrt(np.pi**2*D_rec**2*L_helical_RX/wavelength**3)
        else:
            raise ValueError("wrong anttype input for TX")

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
    pointoff_transm = 0.5
    LF_atten = 0.95
    
    freq = 2.20*10**9
    distance = 1.97*10**6
    
    pointoff_rec = 0.
    D_rec = 15.     
    LF_rec = 0.7
    DataRate = 2.27*10**9
    T_syst = 135.

    SNratio = communication(P_transm,LF_transm,D_transm,pointoff_transm,LF_atten,\
                  freq,distance,pointoff_rec,D_rec,LF_rec,DataRate,T_syst, "dish", "helical",0,1.)
    print SNratio
    
    plt.plot(D_transm,SNratio)
    plt.show()