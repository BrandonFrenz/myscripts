#!/usr/bin/python2.7
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

native = open(sys.argv[1], "r")
natstart = int(sys.argv[2])
natstop = int(sys.argv[3])
newpose = open(sys.argv[4], "r")
posestart = int(sys.argv[5])
posestop = int(sys.argv[6])

natpositions = []
posepositions = []
for line in native.readlines():
    if line.startswith("ATOM"):
       sline = line.split()
       residue = int(sline[5])
       atom = sline[2]
       if(atom == 'CA'):
          if(residue >= natstart and residue <= natstop):
             coords = {};
             x = float(sline[6])
             y = float(sline[7])
             z = float(sline[8])
             coords['x'] = x
             coords['y'] = y
             coords['z'] = z
             natpositions.append(coords)
for line in newpose.readlines():
    if line.startswith("ATOM"):
       sline = line.split()
       residue = int(sline[5])
       atom = sline[2]
       if(atom == 'CA'):
           if(residue >= posestart and residue <= posestop):
              coords = {};
              x = float(sline[6])
              y = float(sline[7])
              z = float(sline[8])
              coords['x'] = x
              coords['y'] = y
              coords['z'] = z
              posepositions.append(coords)
rescount = natstop-natstart+1
distancesum = 0
for i in range(0, rescount):
    nx = natpositions[i]['x']
    ny = natpositions[i]['y']
    nz = natpositions[i]['z']
    px = posepositions[i]['x']
    py = posepositions[i]['y']
    pz = posepositions[i]['z']
    distance = math.sqrt(((nx - px)**2 + (ny - py)**2 + (nz - pz)**2))
    distancesum += distance**2
rms = math.sqrt(distancesum/rescount)
print rms
