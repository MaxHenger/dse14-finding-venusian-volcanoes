# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 14:43:15 2016

@author: MaxHenger
"""

# SomeFunction is called from YuyangExampleMain. You can assume it accepts a lot
# of arguments which are settings. So these settings will not change while the
# program is running.
# This is very important, because if one function modifies the 'configuration'
# variable (which is a class), then it will also be changed outside of it
def SomeFunction(config):
	print('first value =', config.first)
	print('second value =', config.second)

# SomeWrongFunction is called from YuyangExampleMain. It is not wrong per-se,
# but it shows the possible issues that can arise from editing the config file
def SomeWrongFunction(config):
	config.second = "me"
