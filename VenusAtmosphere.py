"""
Created on Mon May 02 17:04:37 2016

@author: Mathijs

Venusian atmospheric model 0-72 km based on 1985 source found in literature study
"""

from scipy import interpolate
from matplotlib import pyplot as plt
import numpy as np

# Enter altitude in km (either one value or an array) and program will return temperature, pressure, density and gravitational acceleration for each altitude
# DO NOT USE ABOVE 72 KM
# CREATED FOR LATITUDE OF 0-30 DEG
# MARGINS STILL TO BE INCLUDED
    
def VenusAtmosphere30latitude(h):
    """Input of altitude in m"""
    h=h/1000.
    # Data for Venus atmosphere with latitude 0-30 deg from 0 to 72 km altitude from 1995 source
    data072 = """0 735.3 92.10 64.79 1.0100 8.06 8.869
    1 727.7 86.45 61.56 1.0083 8.09 8.867
    2 720.2 81.09 58.45 1.0065 8.10 8.864
    3 712.4 76.01 55.47 1.0050 8.12 8.861
    4 704.6 71.20 52.62 1.0034 8.14 8.858
    5 696.8 66.65 49.87 1.0022 8.17 8.855
    6 688.8 62.35 47.24 1.0010 8.19 8.852
    7 681.1 58.28 44.71 1.0000 8.21 8.849
    8 673.6 54.44 42.26 0.9990 8.24 8.846
    9 665.8 50.81 39.95 0.9981 8.26 8.843
    10 658.2 47.39 37.72 0.9973 8.28 8.840
    11 650.6 44.16 35.58 0.9966 8.30 8.837
    12 643.2 41.12 33.54 0.9958 8.32 8.834
    13 635.5 38.26 31.60 0.9954 8.34 8.831
    14 628.1 35.57 29.74 0.9949 8.36 8.829
    15 620.8 33.04 27.95 0.9947 8.39 8.826
    16 613.3 30.66 26.27 0.9944 8.41 8.823
    17 605.2 28.43 24.68 0.9942 8.44 8.820
    18 597.1 26.33 23.18 0.9939 8.46 8.817
    19 589.3 24.36 21.74 0.9937 8.49 8.814
    20 580.7 22.52 20.39 0.9935 8.52 8.811
    21 572.4 20.79 19.11 0.9933 8.55 8.808
    22 564.3 19.17 17.88 0.9931 8.58 8.805
    23 556.0 17.66 16.71 0.9930 8.60 8.802
    24 547.5 16.25 15.62 0.9929 8.64 8.800
    25 539.2 14.93 14.57 0.9929 8.67 8.797
    26 530.7 13.70 13.59 0.9928 8.70 8.794
    27 522.3 12.56 12.65 0.9928 8.74 8.791
    28 513.8 11.49 11.77 0.9928 8.77 8.788
    29 505.6 10.50 10.93 0.9929 8.80 8.785
    30 496.9 9.581 10.15 0.9929 8.84 8.782
    31 488.3 8.729 9.406 0.9930 8.87 8.779
    32 479.9 7.940 8.704 0.9931 8.91 8.776
    33 471.7 7.211 8.041 0.9932 8.95 8.774
    34 463.4 6.537 7.420 0.9933 9.00 8.771
    35 455.5 5.917 6.831 0.9935 9.04 8.768
    36 448.0 5.346 6.274 0.9937 9.08 8.765
    37 439.9 4.822 5.762 0.9939 9.12 8.762
    38 432.5 4.342 5.276 0.9941 9.17 8.759
    39 425.1 3.903 4.823 0.9944 9.22 8.756
    40 417.6 3.501 4.404 0.9946 9.26 8.753
    41 410.0 3.135 4.015 0.9949 9.31 8.750
    42 403.5 2.802 3.646 0.9951 9.35 8.748
    43 397.1 2.499 3.303 0.9954 9.40 8.745
    44 391.2 2.226 2.985 0.9957 9.44 8.742
    45 385.4 1.979 2.693 0.9960 9.48 8.739
    46 379.7 1.756 2.426 0.9962 9.52 8.736
    47 373.1 1.556 2.186 0.9964 9.57 8.733
    48 366.4 1.375 1.967 0.9966 9.62 8.730
    49 358.6 1.213 1.769 0.9988 9.68 8.728
    50 350.5 1.066 1.594 0.9970 9.76 8.725
    51 342.0 0.9347 1.432 0.9971 9.83 8.722
    52 333.3 0.8167 1.284 0.9972 9.91 8.719
    53 323.0 0.7109 1.153 0.9973 10.01 8.716
    54 312.8 0.6160 1.032 0.9974 10.13 8.713
    55 302.3 0.5314 0.9207 0.9975 10.25 8.710
    56 291.8 0.4559 0.8183 0.9976 10.38 8.708
    57 282.5 0.3891 0.7212 0.9978 10.50 8.705
    58 275.2 0.3306 0.6289 0.9979 10.59 8.702
    59 268.7 0.2796 0.5448 0.9981 10.67 8.699
    60 262.8 0.2357 0.4694 0.9982 10.75 8.696
    62 254.5 0.1659 0.3411 0.9986 10.85 8.690
    64 245.4 0.1156 0.2443 0.9989 10.98 8.684
    66 241.0 0.07970 0.1729 0.9993 11.04 8.679
    68 235.4 0.05447 0.1210 0.9994 11.11 8.673
    70 229.8 0.03690 0.08393 0.9995 11.19 8.667
    72 224.1 0.02476 0.05775 0.9996 11.25 8.662""" 
    
    # Split the data
    data  = data072.split()
    # Create empty array
    data_sorted=[]
    # Define number of columns
    n=7
    # Create an array which can contain all 7 sets of data
    for i in range(0,n):
        data_sorted.append([])
    # Add the data to the empty array
    for i,dat in enumerate(data):
        data_sorted[i%n].append(float(dat))
    
    # Define the data
    altlist = data_sorted[0]
    TempK = data_sorted[1]
    PressBar = data_sorted[2]
    Density = data_sorted[3]
    gravAcc = data_sorted[6]
    
    # Interpolate the  data [Checked for accuracy using plots]
    TempK_inter = interpolate.InterpolatedUnivariateSpline(altlist,TempK)
    PressBar_inter = interpolate.InterpolatedUnivariateSpline(altlist,PressBar)
    Density_inter = interpolate.InterpolatedUnivariateSpline(altlist,Density)
    gravAcc_inter = interpolate.InterpolatedUnivariateSpline(altlist,gravAcc)
    
    # Evaluate the values at the designated altitude(s)
    h_Temp = TempK_inter(h) 
    h_PressBar = PressBar_inter(h)
    h_Density = Density_inter(h)
    h_gravAcc = gravAcc_inter(h)
    
    # Return the values
    return (h_Temp,h_PressBar*10**5,h_Density,h_gravAcc)

if __name__=="__main__":
    #altitudes = np.arange(0,73000)
    #Temp, Press, Density, GravAcc = VenusAtmosphere30latitude(altitudes)
    x = np.arange(0,120000)
    y = VenusAtmosphere30latitude(x)
    y2 = scale_height(x)
    
    plt.plot(x,y[2],label="interpolate")
    plt.plot(x,y2[2],label="scale")
    plt.plot(x,np.zeros(len(x)))
    plt.legend()
