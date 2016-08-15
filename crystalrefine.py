#!/usr/bin/python2.7
import sys
from argparse import ArgumentParser
from multiprocessing import Pool
from os import popen, system
from os.path import basename, exists
from sys import exit, stderr, stdout
import random
from pprint import pprint
import operator
import math
import os
import multiprocessing
import time
#takes input as follows:
# i=xml a=pdb b=symfile c=fragfile d=mtz e=nativepdb x=suffix
def runner(i,x):
   # for i in range(7, 12):  
        a = sys.argv[1] 
        b = sys.argv[2]
        c = sys.argv[3]
        d = sys.argv[4]
        e = sys.argv[5]
        h = sys.argv[6]
        os.system("phenix.rosetta.run_phenix_interface ~/Rosetta/main/source/bin/rosetta_scripts.python.linuxgccrelease -database ~/Rosetta/main/database -parser:protocol /work/bfrenz/lowres_model_building/crystal_data/xmlvariants/refine%s.xml -s %s -parser:script_vars symmfile=%s fragfile=%s -cryst::mtzfile %s -crystal_refine -in:file:native %s -out:suffix _%s.%s.%d -nstruct 1 -ignore_unrecognized_res" % (i,a,b,c,d,e,h,i,x) )
        #os.system("sh run.sh refine%s.xml  %s %s %s %s uniqueid_%s.%d" % (i,a,b,c,d,i,x) )

if __name__ == '__main__':         
        jobs = []
        for i in range(7,12):
                for x in range(1, 3):
                        p = multiprocessing.Process(target=runner, args=(i,x,))
                        jobs.append(p)
                        p.start()

           

         
        
