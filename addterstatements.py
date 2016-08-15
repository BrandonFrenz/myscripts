#!/usr/bin/python2.7
import sys
import operator
import math
import os
import time
import glob
import re
import numpy

def main():
    notthesame()

def notthesame():
    pdbs = []
    load_pdbs(pdbs)
    for i in pdbs:
        terminalres = []
        find_gaps(i, terminalres)
        add_ters(i, terminalres)

def allthesame():
    pdbs = []
    terminalres = []
    load_pdbs(pdbs)
    find_gaps(pdbs[0], terminalres)
    print terminalres
    for i in pdbs:
        add_ters(i, terminalres)

def load_pdbs(pdbs):
    it=1
    while it<=len(sys.argv)-1:
        pdbs.append(sys.argv[it])
        it+=1

def find_gaps(name, terminalres):
    with open(name,'r') as firstpdb:
        first = True
        N = numpy.array((0,0,0))
        O = numpy.array((0,0,0))
        resnum = 0
        previousnum = 1
        for line in firstpdb:
            if not line.startswith('ATOM'):
                continue
            line = line.strip()
            line = line.split();
            atom = line[2]
            resnum = int(line[5])
            if atom == 'N' and first == True:
                first = False
            if atom == 'C':
                ox = float(line[6])
                oy = float(line[7])
                oz = float(line[8])
                O = numpy.array((ox, oy, oz))
            if atom == 'N' and first == False:
                nx = float(line[6])
                ny = float(line[7])
                nz = float(line[8])
                N = numpy.array((nx, ny, nz))
            if previousnum != resnum:
                distance = numpy.linalg.norm(N-O)
                if distance > 3.0:
                    terminalres.append(previousnum)
            previousnum = resnum

def add_ters(name, terminalres):
    newpdb = []
    with open(name, 'r') as pdb:
        previousresnum = 1
        resnum = 1
        for line in pdb:
            if line.startswith('ATOM'):
                splitline = line.strip().split()
                resnum = int(splitline[5])
                if resnum != previousresnum and previousresnum in terminalres:
                    newpdb.append('TER\n')
                previousresnum = resnum
            newpdb.append(line)
    with open(name, 'w') as temp:
        for i in newpdb:
            temp.write(i)

main()
