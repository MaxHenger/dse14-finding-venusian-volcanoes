# -*- coding: utf-8 -*-
"""
Created on Tue May 31 12:11:35 2016

@author: Chaggai
"""
import os.path
files = [
         
         ]
 
checked = {}
 
def analyze(name):
    name = name.strip()
    if os.path.exists(name+".py"):
        if name+".py" not in files:
            files.append(name+".py")
 
i = 0
n = 0
while i<len(files):
    name = files[i]
#     print name
    checked[name] = 0
    f = open(name)
    docOpened = False
    for j,line in enumerate(f):
        if not docOpened and "'''" in line:
            if line.count("'''")==1:
                docOpened = True
            continue
        elif docOpened:
            if "'''" in line:
                docOpened = False
            continue
       
        if line.lstrip().startswith("#") or line.strip()=="":
            continue
        n+=1
       
        if "import" in line:
            line = line.lstrip("import")
            line = line.lstrip(" ")
            if "," in line:
                names = map(str.strip,line.split(","))
                for name in names:
                    analyze(name)
            else:
                analyze(line)
    i+=1
               
print "Found {} lines of code".format(n)