# -*- coding: utf-8 -*-
"""
Created on Fri May 27 14:29:46 2016

@author: Chaggai
"""
#lander-orbiter = 40.7 hr/7.6hr

#orbiter-aircraft = 
#orbiter-earth = 

def datastorage_AC(dataratepayload, dataratesystem, contacttime, P_SC, frac_sciencetime, safetyfactor, mean_compres_fac=1.):
    #datarates in bits/s
    t_nocontact = P_SC-contacttime
    t_science = t_nocontact*frac_sciencetime
    Datahousekeeping = dataratesystem*P_SC/8. #bytes
    Datapayload = dataratepayload*t_science/8. #bytes
    totalData= safetyfactor*((Datahousekeeping+Datapayload)/mean_compres_fac) #bytes
    totalDatagig = totalData/(1.*10**9) #total data in GByte
    
    #mass solid state recorder
    M_recorder = totalDatagig/(0.041*totalDatagig+0.3128) #relation taken from zandbergen lecture notes AE1222 (outdated just for inital use)
    
    return(totalData, M_recorder)
    
def datastorage_L(dataratepayload, dataratesystem, lifetime_L, frac_sciencetime, safetyfactor, mean_compres_fac=1.):
    #datarates in bits/s
    t_science = lifetime_L*frac_sciencetime
    Datahousekeeping = dataratesystem*lifetime_L/8. #bytes
    Datapayload = dataratepayload*t_science/8. #bytes
    totalData= safetyfactor*((Datahousekeeping+Datapayload)/mean_compres_fac) #bytes
    totalDatagig = totalData/(1.*10**9) #total data in GByte
    
    #mass solid state recorder
    M_recorder = totalDatagig/(0.041*totalDatagig+0.3128) #relation taken from zandbergen lecture notes AE1222 (outdated just for inital use)
    
    return(totalData, M_recorder)
    
def datastorage_SC(dataratepayloadSC,dataratesystemSC, DataAC, DataL, contacttime_earth, P_SC, frac_sciencetime, safetyfactor, mean_compres_fac=1.): 
    #datarates in bits/s
    t_nocontact_earth = P_SC-contacttime_earth
    t_science = t_nocontact_earth*frac_sciencetime
    DatahousekeepingSC = dataratesystemSC*P_SC/8. #bytes
    DatapayloadSC = dataratepayloadSC*t_science/8. #bytes
    totalDataSC= safetyfactor*((DatahousekeepingSC+DatapayloadSC)/mean_compres_fac)+DataAC+DataL #bytes
    totalDatagigSC = totalDataSC/(1.*10**9) #total data in GByte
    
    #mass solid state recorder
    M_recorderSC = totalDatagigSC/(0.041*totalDatagigSC+0.3128) #relation taken from zandbergen lecture notes AE1222 (outdated just for inital use)
    
    return(totalDataSC, M_recorderSC)
    
    
if __name__ == "__main__":
    
    import communications as com
    import numpy as np
    secundairy =True    
    
    safetyfactor = 1.5
    mean_compres_fac = 15.0    
    P_SC = 2894.9*60.
    
    #aircraft primary ins
    camdatarate1_AC = 3840.*2160.*24/60.
    camdatarate2_AC = 1080.*1920.*24/120.
    VEMdatarate_AC = 190.*1000.
    householddatarate_AC = 600.*8    
    #secondary aircraft
    spectro_AC = 1.46*1000.
    spectro_xray_AC = 32*1000.
    if secundairy ==True:
        dataratepayload_AC = camdatarate1_AC+camdatarate2_AC+VEMdatarate_AC+spectro_AC+spectro_xray_AC
    else:
        dataratepayload_AC = camdatarate1_AC+camdatarate2_AC+VEMdatarate_AC    
    contacttime_AC = 40.7*3600.
    frac_sciencetime_AC = 0.6
    dataratebud_AC = (dataratepayload_AC*frac_sciencetime_AC+householddatarate_AC)*P_SC/contacttime_AC/mean_compres_fac*safetyfactor   #for the link budget
    
    #Lander
    Chemcam = 1032.*1000./30.
    camdatarate1_L = 3840.*2160.*24/100.
    IRdatarate_L = 1080.*1920.*12/30.
    householddatarate_L = 600.*8.
    dataratepayload_L = camdatarate1_L+IRdatarate_L+householddatarate_L
    lifetime_L = 3.*3600.
    frac_sciencetime_L = 1.
    dataratebud_L = (dataratepayload_L+householddatarate_L)/mean_compres_fac*safetyfactor #for the link budget


    #spacecraft
    SARdatarate_SC = 800.*1000.
    VEMdatarate_SC = 190.*1000.
    householddatarate_SC = 800.*8.    
    #spacecraft secondary
    Lidar_SC = 80.*1000. #assumption   
    camera_SC=1920*1080*24/1200. 
    if secundairy ==True:
        dataratepayload_SC = SARdatarate_SC+VEMdatarate_SC+householddatarate_SC+Lidar_SC+camera_SC
    else: 
        dataratepayload_SC = SARdatarate_SC+VEMdatarate_SC+householddatarate_SC
    contacttime_SC = 7.6*3600.
    frac_sciencetime_SC = 1.
    dataratebud_SC = (dataratepayload_SC*frac_sciencetime_SC+householddatarate_SC)*P_SC/contacttime_SC/mean_compres_fac*safetyfactor+dataratebud_L/6.+dataratebud_AC #for link budget
    
    
    Data_AC, M_recorder_AC =  datastorage_AC(dataratepayload_AC, householddatarate_AC, contacttime_AC, P_SC, frac_sciencetime_AC, safetyfactor, mean_compres_fac)
    Data_L, M_recorder_L =  datastorage_L(dataratepayload_L, householddatarate_L, lifetime_L, frac_sciencetime_L, safetyfactor, mean_compres_fac)
    Data_SC, M_recorder_SC = datastorage_SC(dataratepayload_SC,householddatarate_SC, Data_AC, Data_L, contacttime_SC, P_SC, frac_sciencetime_SC, safetyfactor, mean_compres_fac)
    
    print "Aircraft", Data_AC/(1.*10**9), M_recorder_AC
    print "Lander", Data_L/(1.*10**9), M_recorder_L
    print "Spacecraft", Data_SC/(1.*10**9), M_recorder_SC
    
    
    
    
    
    
    
    #LINK BUDGETS
    turnaround = {"sband": 221./240., "xband": 749./880., "kaband": .93}
    bands={"sband":2.23*10**9, "xband":8.4*10**9,"kaband":34.5*10**9}
    tsys = {"sband": 135., "xband": 135., "kaband": 424.}
    coding = [5., 4., 7./2., 3, 8./3., 5./2., 7./3., 9./4., 11./5., 17./8., 19./9., 8./3.,5./2.,7./3.,11./5.,17./8.,19./9.,5./2.,7./3.,9./4.,11./5.,17./8.,19./9.]
    codingname = ["QPSK1/4", "QPSK1/3","QPSK2/5","QPSK1/2","QPSK3/5","QPSK2/3","QPSK3/4","QPSK4/5","QPSK5/6","QPSK8/9","QPSK9/10", "8PSK3/5","8PSK2/3","8PSK3/4","8PSK5/6","8PSK8/9","8PSK9/10",
                  "16APSK2/3","16APSK3/4","16APSK4/5","16APSK5/6","16APSK8/9","16APSK9/10"]
    codemargin = [0.75,.59,.73,1.05,1.48,1.89,2.31,2.67,2.99,3.73,3.89,3.0,3.65,4.43,5.41,6.46,6.7,4.76,5.49,6.03,6.42,7.42,7.61]
    
    
    #spacecraft parameters to AC and lander
    Ptrans_SC_L_AC = 20. #watts
    LF_SC_L_AC = 0.8
    D_SC_L_AC = 0.8 #m
    pointoff_SC_L_AC = 0.2 #degrees
    anttype_SC_L_AC = "horn"
    
    #groundstation charecteristics
    Ptrans_Ground = 20000. #watts
    LF_Ground = 0.8
    D_Ground = 70. #m
    pointoff_Ground = 0.006 #degrees http://ipnpr.jpl.nasa.gov/progress_report/42-87/87W.PDF
    anttype_Ground = "horn"
    
    #lander-SC 
    band_L = "xband"
    turnaround_L = turnaround[band_L]
    landdlink_L = bands[band_L]
    landulink_L = turnaround_L*landdlink_L
    T_sys_L = tsys[band_L]
    Ptrans_L = 28. #watts    
    LF_L = 0.8    
    anttype_L = "horn"    
    D_L = 0.26#m    
    pointoff_L = 0.5 #degrees    
    LF_atten_L = 1.
    dist_L = 112022.5*1000. #m
    
    #aircraft-SC 
    band_AC = "xband"
    turnaround_AC = turnaround[band_AC]
    landdlink_AC = bands[band_AC]
    landulink_AC = turnaround_AC*landdlink_AC
    T_sys_AC = tsys[band_AC]
    Ptrans_AC = 23. #watts
    LF_AC = 0.8
    anttype_AC = "dish"
    D_AC = 0.28#m
    pointoff_AC = 0.5 #degrees
    LF_atten_AC = 1.
    dist_AC = 111992.5*1000. #m
    
    #SC -Earth
    band_SC = "kaband"
    turnaround_SC = turnaround[band_SC]
    landdlink_SC = bands[band_SC]
    landulink_SC = turnaround_SC*landdlink_SC
    T_sys_SC = tsys[band_SC]
    Ptrans_SC = 280. #watts
    LF_SC = 0.8
    anttype_SC = "dish"
    D_SC = 2.5#m
    pointoff_SC = 0.006 #degrees
    LF_atten_SC = 0.95
    dist_SC = 275810000.*1000. #m    
    
    
    
    Downlink_L=[]
    for count in range(len(coding)):
        SNRdown_L = com.communication(Ptrans_L,LF_L,D_L,pointoff_L,LF_atten_L,\
                  landdlink_L,dist_L,pointoff_SC_L_AC,D_SC_L_AC,LF_SC_L_AC,coding[count]*dataratebud_L,T_sys_L, anttype_L, anttype_SC_L_AC,L_helical_TX=0, L_helical_RX=0)
        Downlink_L.append(SNRdown_L-6.-codemargin[count])
    index_L=np.argmax(Downlink_L)
    print "lander downlink", codingname[index_L],Downlink_L[index_L]     
    
    
    Downlink_AC=[]
    for count in range(len(coding)):
        SNRdown_AC = com.communication(Ptrans_AC,LF_AC,D_AC,pointoff_AC,LF_atten_AC,\
                  landdlink_AC,dist_AC,pointoff_SC_L_AC,D_SC_L_AC,LF_SC_L_AC,coding[count]*dataratebud_AC,T_sys_AC, anttype_AC, anttype_SC_L_AC,L_helical_TX=0, L_helical_RX=0)
        Downlink_AC.append(SNRdown_AC-6.-codemargin[count])
    index_AC=np.argmax(Downlink_AC)
    print "aircraft downlink", codingname[index_AC],Downlink_AC[index_AC]    
    
    Downlink_SC=[]
    for count in range(len(coding)):
        SNRdown_SC = com.communication(Ptrans_SC,LF_SC,D_SC,pointoff_SC,LF_atten_SC,\
                  landdlink_SC,dist_SC,pointoff_Ground,D_Ground,LF_Ground,coding[count]*dataratebud_SC,T_sys_SC, anttype_SC, anttype_Ground,L_helical_TX=0, L_helical_RX=0)
        Downlink_SC.append(SNRdown_SC-3.-codemargin[count])
    index_SC=np.argmax(Downlink_SC)
    print "spacecraft downlink", codingname[index_SC],Downlink_SC[index_SC]   
    
    
    
    