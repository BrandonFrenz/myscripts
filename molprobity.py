#!/usr/local/bin/python2.7
import argparse
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument('-s','--structures', nargs='+', help="the list of structures to be renamed")
args = parser.parse_args()
pdbs = args.structures
for i in pdbs:
    molprobity = 'phenix.molprobity %s'%i
    os.system(molprobity)
    name = re.sub('.pdb','',i)+'mpscore'
    os.system('mv molprobity.out %s'%name)
