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


argone = sys.argv[1] #xml
argtwo = sys.argv[2] #input pdb
argthree = sys.argv[3] #symm
argfour = sys.argv[4] #path to fragfile
argfive = sys.argv[5] #mtz file
argsix = sys.argv[6] #native pdb

#make density map
os.system("~bfrenz/scripts/make_map_phenix.sh %s %s 1.0 " % (argtwo,argfive) )

#make file with density correlation 
os.system("/work/bfrenz/Rosetta/main/source/bin/density_tools.default.linuxgccrelease -s %s -mapfile inputdensity.map -perres -mute all &> densitycorrelation.txt" % (argtwo) )

#get density correlation and residue list
file = open("densitycorrelation.txt", "r")

cor_list = []; res_list = [];
for line in file.readlines():
      if line.startswith("residue"):
          l = line.split()
          cor = float( l[2].split("=")[1] )
          resnum = int( l[1] )
          cor_list.append(cor)
          res_list.append(resnum)
file.close()                
zscore = stats.zscore(cor_list)
#get targets
def movingaverage(interval, window_size):
    window = numpy.repeat(1.0, window_size)/(window_size)
    smas = numpy.convolve(interval,window,'valid')
    return smas #
windowsize = 20 #change to vary smoothing.
corav = movingaverage((zscore),windowsize) 
minima = []
coraver = corav.tolist()
loopnum = 3 #how many loops do you want to rebuild?
for z in range(0, windowsize):
      coraver.append(99)
      coraver.insert(0, 99)
for i in range(0, loopnum):
   ca = coraver.index(min(coraver))
   if ca < windowsize+19:
       minima.append(0)
       minima.append(19)
   else:
       min20 = res_list.index(ca-20) 
       min0 = res_list.index(ca)
       minima.append(min20)
       minima.append(min0)
   coraver[ca] = '999' 
   for x in range(0, windowsize):      
         y = ca+x
         z = ca-x
         coraver[y] = '9999'   
         coraver[z] = '9999'
         coraver[ca] = '9999'
file = open(argtwo)
for line in file.readlines():
    if line.startswith("ATOM"):
        l = line.split()
        chainname = str( l[4] )
for i in range(0,6,2):         
      os.system(
        "phenix.rosetta.run_phenix_interface \
            ~/Rosetta/main/source/bin/rosetta_scripts.python.linuxgccrelease \
            -database ~/Rosetta/main/database \
            -parser:protocol %s \
            -s %s \
            -parser:script_vars symmfile=%s fragfile=%s pos1=%s%s pos2=%s%s \
            -cryst::mtzfile %s \
            -crystal_refine \
            -in:file:native %s \
            -out:suffix _%s  \
            -nstruct 1 \
            -ignore_unrecognized_res" % (argone,argtwo,argthree,argfour,minima[i],chainname,minima[i+1],chainname,argfive,argsix,i) )

      

