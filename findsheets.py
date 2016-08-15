#!/usr/local/bin/python2.7
import sys
from argparse import ArgumentParser
from multiprocessing import Pool
from os import popen,  system
from os.path import basename, exists
from sys import exit, stderr, stdout
import operator
import os
import multiprocessing
import time
import glob

beams = glob.glob("beam*")
for files in beams:
    hassheets = False
    beam = open("%s" % files, "r")
    for line in beam:
        splitline = line.split()
        if len(splitline) == 3:
            hassheets = True
    if hassheets:
        print files
