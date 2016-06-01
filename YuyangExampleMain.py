# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 14:42:52 2016

@author: MaxHenger
"""

import YuyangExampleOther as yeo

# Config is a class that contains all kinds of variables. You can see a class
# as a type of container: it is a variable, but it can hold variables itself.
# There are two ways of setting variables. The way I personally like is the
# one that is being used, the second one is commented out.
# The reason why I like to use the first one (in the __init__ function) is also
# shown in the main progam execution
class Config:
	#first = "hello"
	#second = "world"

	def __init__(self, setupConfigVariable=0):
		# The keyword 'self' indicates that you're modifying variables belonging
		# to this class.
		self.first = "hello"
		self.second = "world"

		if setupConfigVariable == 1:
			self.second = "universe"

# The main script starts here
# - first we create an instance of the 'Config' class. You can have multiple
# 	instances lying around. Note that each instance of the config class is its
# 	own container.
conf = Config() # this will create an instance of Config. It will call Config's __init__ function

# - We call the function that requires a large set of configuration variables
print(' * Simply calling a function accepting the Config class')
yeo.SomeFunction(conf)

# - Now we create a second instance of the Config class and call the same
# 	same function. However, now we use the auxilliary variable to change the
#	way in which the class is initialized
print(' * Creating an alternative instance of the Config class')
confAlternative = Config(1) # see how setupConfigVariable is now set to 1

# - We call the function first on the alternative configuration, then on the
#	original configuration. See how both 'containers' hold different instances
#	of their member variables?
yeo.SomeFunction(confAlternative)
yeo.SomeFunction(conf)

# - Now we show the dangers of editing values in the Config class in a function
#	while expecting it to stay exactly the same. We will change the first
#	instance of the config class
print(' * Changing the config class from within a function')
yeo.SomeWrongFunction(confAlternative)
yeo.SomeFunction(confAlternative)
yeo.SomeFunction(conf)

# So, what is useful about the config class? Lets say you need to add a shitload
# of parameters to a large number of functions, that will require a lot of
# typing! If you instead define such a helper-class, as defined above, you only
# have to pass in one argument each time. Not only that: If you decide to change
# all of your functions to depend on different variables, it will only require
# one change to the helper-class.
