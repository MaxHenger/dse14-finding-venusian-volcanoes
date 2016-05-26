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
    def __init__(self,definition="complex",accuracy=10,uncertainty=False):
        if definition=="complex":
            self.accuracy=accuracy
            self.uncertainty=uncertainty
            self.__importJ__(self.accuracy,self.uncertainty)
        else:
            raise ValueError("Unknown Definition")
    
    def __importJ__(self,accuracy,uncertainty,filePath="./shgj180p.a01"):
        self.accuracy=accuracy
        self.uncertainty=uncertainty
        with open(filePath,"r") as const:
            lines = const.readlines()
            
            header=lines[237-1]
            header=header.split(",")
            headerLi=[]
            for i in header:
                headerLi.append(i.strip())
            headerLi[0]=headerLi[0].split()[1]
            #headerMeaning=["Name","Radius","GM","unc GM","degree","order","Normalization state","Unknown","Unknown"]
            if uncertainty:
                self.Mu = np.array([float(headerLi[1])-float(headerLi[2]),float(headerLi[1]),float(headerLi[1])+float(headerLi[2])])*10**9
                self.RadiusMean = np.array([float(headerLi[0])]*3)*1000
            else:
                self.Mu = np.array(float(headerLi[1]))*10**9
                self.RadiusMean = np.array(float(headerLi[0]))*1000
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
                self.C[int(cofloat[0])][int(cofloat[1])]=float(cofloat[2])
                self.S[int(cofloat[0])][int(cofloat[1])]=float(cofloat[3])
                self.uncC[int(cofloat[0])][int(cofloat[1])]=float(cofloat[4])
                self.uncS[int(cofloat[0])][int(cofloat[1])]=float(cofloat[5])
            self.C[0][0]=1.0
            
            
class GravityEarth:
    def __init__(self,definition="complex",accuracy=10):
        if definition=="complex":
            self.accuracy=accuracy
            self.__importJ__(self.accuracy)
        else:
            raise ValueError("Unknown Definition")
            
    def __importJ__(self,accuracy,filePath="./egm96.a01"):
        self.RadiusMean=6378137.0 
        self.Mu=0.3986004418e+15
        self.accuracy=accuracy
        with open(filePath,"r") as const:
            lines = const.readlines()
            if accuracy==0:
                size = (360+2,360+2)
            else:
                size = (accuracy,accuracy)
            self.C=np.zeros(size)
            self.S=np.zeros(size)
            self.uncC=np.zeros(size)
            self.uncS=np.zeros(size)
            coeff = lines
            for co in coeff:
                co=co.split()
                cofloat=[]
                for i,value in enumerate(co):
                    cofloat.append((value.strip()))
                if float(cofloat[0])>accuracy-1:
                    break
                self.C[int(cofloat[0])][int(cofloat[1])]=float(cofloat[2])
                self.S[int(cofloat[0])][int(cofloat[1])]=float(cofloat[3])
                self.uncC[int(cofloat[0])][int(cofloat[1])]=float(cofloat[4])
                self.uncS[int(cofloat[0])][int(cofloat[1])]=float(cofloat[5])
    
class GravitySun:
    def __init__(self):
        pass