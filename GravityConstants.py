# -*- coding: utf-8 -*-
"""
Created on Tue May 10 19:46:30 2016
This file contains the constants that are used inside the gravity modeling
code in Gravity.py (defining the Gravity class). This file does not need any 
editing if one is not attempting to edit or add a new model.
I know this file is very ... sparse, however if we want to develop this further
it will serve as a good basis.
@author: Julius
"""
import numpy as np

class GravityVenus:
    def __init__(self,definition="Magellen",accuracy=10):
        if definition=="simple":
            self.Mu         = np.array([0.32486 *10**6 *10**9]*3)
            self.RadiusMean = np.array([6051.8 *1000]*3)
            self.J          = np.zeros((3,3))
            self.uncJ       = np.zeros((3,3))
            self.lam        = np.zeros((3,3))
            self.unclam     = np.zeros((3,3))
            self.J[2][0]    = 4.458 *10**-4  # the J2 effect
        elif definition=="Magellen":
            self.accuracy=accuracy
            self.__importJ__(self.accuracy)
        else:
            raise ValueError("Unknown Definition")
    
    def __importJ__(self,accuracy,filePath="./shgj120p.a01"):
        self.accuracy=accuracy
        with open(filePath,"r") as const:
            lines = const.readlines()
            
            header=lines[237-1]
            header=header.split(",")
            headerLi=[]
            for i in header:
                headerLi.append(i.strip())
            headerLi[0]=headerLi[0].split()[1]
            headerMeaning=["Name","Radius","GM","unc GM","degree","order","Normalization state","Unknown","Unknown"]
            self.Mu = np.array([float(headerLi[1])-float(headerLi[2]),float(headerLi[1]),float(headerLi[1])+float(headerLi[2])])*10**9
            self.RadiusMean = np.array([float(headerLi[0])]*3)*1000
            coeff = lines[238-1:]
            if accuracy==0:
                size = (int(headerLi[3])+2,int(headerLi[4])+2)
            else:
                size = (accuracy,accuracy)
                
            self.C=np.zeros(size)
            self.S=np.zeros(size)
            self.uncC=np.zeros(size)
            self.uncS=np.zeros(size)
            for co in coeff:
                co=co.split(",")
                cofloat=[]
                for i,value in enumerate(co):
                    cofloat.append((value.strip()))
                if float(cofloat[0])>accuracy-1:
                    break
                self.C[float(cofloat[0])][float(cofloat[1])]=float(cofloat[2])
                self.S[float(cofloat[0])][float(cofloat[1])]=float(cofloat[3])
                self.uncC[float(cofloat[0])][float(cofloat[1])]=float(cofloat[4])
                self.uncS[float(cofloat[0])][float(cofloat[1])]=float(cofloat[5])
                
class GravityEarth:
    def __init__(self):
        self.J          = np.zeros((7,7))
        self.lam        = np.zeros(self.J.size)
        self.RadiusMean = 6378*1000.
        self.Mu         = 398600.4415*10**9
        self.J[2][0]    = 1082.6267 *10**-6 # the J2 effect
        self.J[3][0]    = -2.5327 *10**-6
        self.J[4][0]    = -1.6196 *10**-6 
        self.J[5][0]    = -0.2273 *10**-6 
    
class GravitySun:
    def __init__(self):
        pass
Status API Training Shop Blog About
© 2016 GitHub, Inc. Terms Privacy Security Contact Help