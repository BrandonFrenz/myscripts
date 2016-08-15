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

beams = sys.argv[1]
beams = open("%s" % (beams),"r")
for line in beams.readlines():
    if len(line.split()) == 2:
       newline = line.split()
       print newline[0]

