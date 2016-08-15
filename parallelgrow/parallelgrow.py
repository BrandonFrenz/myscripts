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

xml = sys.argv[1] #xml
inputpdb = sys.argv[2] #input pdb
densitymap = sys.argv[3] #map
native = sys.argv[4] #Native
compnumber = int(sys.argv[5]) # numberof processors

counter = 0
print len(sys.argv)
if len(sys.argv) >6:
    counter = int(sys.argv[6])
filecount = len(glob.glob("beam_%s_*.txt" % (counter),))
if counter == 0:
    counter = 1
beamcount = 0
if filecount == 0:
    filecount = 1
print filecount

def run(filename, jobcounter):
     os.system("sh prun.sh %s %s %s %s 1 %s 1 %s" % (xml,inputpdb,densitymap,native,filename,jobcounter) )
     print " just ran on " + filename

def combine():
    combinedbeams = open("beam_%d.txt" % (counter), "w")
    for i in range(1, filecount+1):
        singleroundbeam = open("beam_%d.%d.txt" % (counter, i), "r")
        combinedbeams.write(singleroundbeam.read())
        singleroundbeam.close()
    combinedbeams.close()
 
def secondfilter(): 
    global beamcount
    filename = ("beam_%d.txt" % (counter))
    os.system("sh prun.sh %s %s %s %s 1 %s 1337 0" % (xml,inputpdb,densitymap,native,filename) )
    print "finished sorting " + filename
    storedbeam = open("beam_0.txt", "r")
    beamcount = 0
    for line in storedbeam.readlines():
        if len(line.split()) == 5: 
            beamcount = beamcount + 1
 
    storedbeam.close()

def splitjobs():
    global filecount
    storedbeam = open("beam_0.txt", "r")
    divisioncount = 1
    partialset = open("beam_0.txt", "r")
    for line in storedbeam.readlines():
        if len(line.split()) == 5:
            partialset.close()
            partialset = open("beam_%d_%d.txt" % (counter,divisioncount), "w")
            divisioncount = divisioncount + 1
        if divisioncount > compnumber:
            divisioncount = 1
        partialset.write(line)
    storedbeam.close()    
    partialset.close()   
    filecount = beamcount
    if beamcount > compnumber :
         filecount = compnumber

def runparalleljobs():
    jobcount = 1
    global filecount
    while jobcount <= filecount:
          if __name__ == '__main__':
             jobs = []
             for i in range(1, compnumber+1):
                 p = multiprocessing.Process(target=run, args=("beam_%d_%d.txt" % (counter,i),i))
                 jobs.append(p)
                 p.start() 
                 jobcount = jobcount + 1
                 if jobcount > filecount:
                     break
             for x in jobs:
                 x.join()


start_time = time.time()
os.system("sh prun.sh %s %s %s %s 0 NA 1 1" % (xml,inputpdb,densitymap,native) )
print "after first"
done = False
while done == False:
    combine()
    secondfilter()
    splitjobs()
    runparalleljobs()
    counter = counter + 1
    if(os.path.isfile("./finished.txt")):
       done = True
       print "found finished"

totaltime = time.time()-start_time
print "the total time to run is "
print totaltime
