# -*- coding: utf-8 -*-
"""
Created on Tue May 17 16:57:23 2016

@author: MaxHenger
"""

import AerodynamicsWing as aerowing
import AerodynamicsUtil as aeroutil

class Analysis:
    def AddWing(self, wing):
        raise RuntimeError("AddWing() is called in base class 'Analysis'")