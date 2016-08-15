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

##This script will calculate the RMSD of all the fragments in a fragfile compared to the torsions given
#to it in a file containing the residue positioni, phi, psi, and omega, in that order
#Enter the fragfile as the first argument (either fragfile format works) and the torsion file as the second argument
positions = {'fragnumbers': {'fragresidue' : { 'torsions'}}};
fragfile = open(sys.argv[1], "r")
nextfrag = True;
fragresidue = {};  fragmentsize = 0; fragnumber = 0;
fragcount = [];
##read in fragment torsions and store in dictionary
i = 0
with fragfile:
        first_line = fragfile.readline()
        if first_line.startswith( "FRAME" ):
          filetype = "newformat"
        elif first_line.startswith( "position" ):
          filetype = "oldformat"
        else:
          filetype = "wrongformat"
fragfile = open(sys.argv[1], "r")
if filetype == "newformat":
    print "new frag file format"
    firstline = True
    for line in fragfile.readlines():   
        a = False
        b= False
        l = line.split()
        if line.startswith( "FRAME" ):
              resnumber = int ( l[1] )
              if(firstline == False):
                  fragcount.append(fragnumber)
              firstline = False    
              fragnumber = 1
              a = True
              fragnumbers = {};
        if line in ['\n', '\r\n']:
            if ( nextfrag == True ):
              fragnumber = fragnumber + 1
              nextfrag = False
            b = True
            i=0
            fragresidue = {};
        if (a == False and b == False):
            nextfrag = True
            i = i + 1
            phi = float(l[6])
            psi = float(l[7])
            omega = float(l[8])
            fragmentsize = i-1
            torsions = {}; 
            torsions['phi'] = phi
            torsions['psi'] = psi
            torsions['omega'] = omega
            fragresidue[i]= torsions
            fragnumbers[fragnumber] = fragresidue
            positions[resnumber] = fragnumbers
if filetype == "oldformat":
  print "old frag file format"  
  firstline = True
  for line in fragfile.readlines():
    a = False
    b= False
    l = line.split()
    if line.startswith( "position" ):
          if(firstline == False):
              fragcount.append(fragnumber)
          firstline = False
          resnumber = int ( l[1] )
          fragnumber = 0
          a = True
          fragnumbers = {};
    if line in ['\n', '\r\n']:
          fragnumber = fragnumber + 1
          b = True
          i=0
          fragresidue = {};
    if (a == False and b == False):
        i = i + 1
        phi = float(l[5])
        psi = float(l[6])
        omega = float(l[7])
        fragmentsize = i-1
        torsions = {}; 
        torsions['phi'] = phi
        torsions['psi'] = psi
        torsions['omega'] = omega
        fragresidue[i]= torsions
        fragnumbers[fragnumber] = fragresidue
        positions[resnumber] = fragnumbers
fragcount.append(fragnumber)        
if filetype == "wrongformat":
    print "##### WRONG FRAG FILE FORMAT ####"
    quit() 
#pprint(positions)
#read in native torsions
ntors = {}; nres = {}; residuenumbering = [];        
nativetorsions = open(sys.argv[2], "r")
for line in nativetorsions.readlines():
    k = line.split()
    nresnum = int(k[0])
    residuenumbering.append(nresnum)
    nphi = float(k[1])
    npsi = float(k[2])
    nomega = float(k[3])
    ntors = {}
    ntors['phi'] = nphi
    ntors['psi'] = npsi
    ntors['omega'] = nomega
    nres[nresnum] = ntors
#pprint(nres)
print "finished read"
#calculate rms
plot = open("fragtorsionrms.txt", "w")
for i in range (1, int (nresnum) - int (fragmentsize)): #for every position
 fragnumber = fragcount[i-1];
 print fragnumber
 for j in range(1, int (fragnumber-1)): #for every fragment at that position 
   rmsd = 0; rmsdnum = 0;   
   for x in range (1, fragmentsize+1): #for every residue in that fragment
    nphi = nres[i+x-1]['phi']
    npsi = nres[i+x-1]['psi']
    nomega = nres[i+x-1]['omega']
    fphi = positions[i][j][x]['phi']      
    fpsi = positions[i][j][x]['psi']
    fomega = positions[i][j][x]['omega']
    if (nphi <0):
      nphi = 180+nphi
    if (npsi <0):
      npsi = 180+npsi
    if (nomega <0):
      nomega = 180+nomega
    if (fphi <0):
      fphi = 180+fphi 
    if (fpsi <0):
      fpsi = 180+fpsi
    if (fomega <0):
      fomega = 180+fomega
    
    rmsdnum = rmsdnum + ((nphi - fphi)**2 + (npsi-fpsi)**2 + (nomega - fomega)**2)
    #rmsd = math.sqrt((((nres[i]['phi'] - (positions[i][j][x]['phi']))**2)))
    if (x == fragmentsize):
         rmsd = math.sqrt(rmsdnum/(fragmentsize*3))  
         numbering = str(i)
         rootmean = str(rmsd)
         plot.write(numbering)
         plot.write(" ")
         plot.write(rootmean)
         plot.write("\n")
