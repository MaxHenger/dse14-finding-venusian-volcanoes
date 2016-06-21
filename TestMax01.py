#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 01:42:13 2016

@author: MaxHenger
"""

class SuperClass:
    def __init__(self, a):
        self.a = a
        
    def method(self, arg1, arg2):
        return arg1 * self.a + arg2
        
class SubClassImplements1(SuperClass):
    def __init__(self, a):
        super().__init__(a)
    
    def method(self, arg1, *args, **kwargs):
        print('*args:', args)
        print('*kwargs:', kwargs)
        return arg1 * self.a
        
class SubClassImplements2(SuperClass):
    def __init__(self, a):
        super().__init__(a)
        
    def method(self, arg1, arg2):
        return arg1 * self.a - arg2
        
cS = SuperClass(5.0)
c1 = SubClassImplements1(5.0)
c2 = SubClassImplements2(5.0)

print(cS.method(2.0, 3.0))
print(c1.method(2.0, 3.0))
print(c2.method(2.0, 3.0))