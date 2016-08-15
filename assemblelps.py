#!/usr/local/bin/python2.7
import sys
from multiprocessing import Pool
from os import popen,  system
from os.path import basename, exists
from sys import exit, stderr, stdout
import random
from pprint import pprint
import operator
import os
import multiprocessing
import time
import glob

loopcount = len(sys.argv)-1 #number of loops
it = 1
lpsfiles = []
while(it < len(sys.argv)):
    lpsfiles.append(sys.argv[it])
    it+=1

combinedlps = open("lpsfile.txt", "w")
stringloopcount = str(loopcount)
combinedlps.write(stringloopcount+ "\n")
for i in lpsfiles:
    beamcount = 0;
    lpsfile = open(i,"r")
    for line in lpsfile.readlines():
        if len(line.split()) == 1:
            beamcount = beamcount + 1
    beamcount = beamcount/2        
    strbeamcount = str(beamcount)
    combinedlps.write(strbeamcount + "\n")
    lpsfile.seek(0)
    combinedlps.write(lpsfile.read())
 
