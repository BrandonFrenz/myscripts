#!/usr/local/bin/python2.7
import sys
from argparse import ArgumentParser
from multiprocessing import Pool
from os import popen,  system
from os.path import basename, exists
from sys import exit, stderr, stdout
import random
from pprint import pprint
import operator
import math
import os
import multiprocessing
import time
import numpy
from scipy import stats
import glob

loopcount = int(sys.argv[1]) #number of loops
combinedlps = open("lpsfile.txt", "w")
stringloopcount = str(loopcount)
combinedlps.write(stringloopcount+ "\n")
for i in range(1, loopcount+1):
    beamcount = 0;
    for x in glob.glob('%s/lpsfile*' % (i)):
        lpsfile = open(x,"r")
        for line in lpsfile.readlines():
            if len(line.split()) == 5:
                beamcount = beamcount + 1
        strbeamcount = str(beamcount)
        combinedlps.write(strbeamcount + "\n")
        lpsfile.seek(0)
        combinedlps.write(lpsfile.read())

