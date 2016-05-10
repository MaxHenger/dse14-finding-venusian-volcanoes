# -*- coding: utf-8 -*-
"""
Created on Tue May 10 09:51:06 2016

@author: Yuyang
"""

import matplotlib.pyplot as plt
from numpy import *

p = arange(300*10**5, 750*10**5, 10**5)
cf = 1.1
r = 4124.
m = 71.

vlist = []
t = 293

rho = p/cf/r/t
v = m/rho

print rho[399]
print v[399]

plt.plot(p, v)
plt.show()