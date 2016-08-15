#!/usr/local/bin/python2.7
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-n','--name', help="the new basename of the files")
parser.add_argument('-s','--structures', nargs='+', help="the list of structures to be renamed")
args = parser.parse_args()
name = args.name
pdbs = args.structures
count = 0
for i in pdbs:
    count+=1
    print i
    os.system(' mv %s %s%s.pdb ' % (i,name,count) )
