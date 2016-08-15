#!/usr/bin/python
import argparse
import pdbtools
import re
import operator
import os

def main():
    args = parseargs()
    molmap = get_mp(args)
    select_top(args,molmap)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdbs',nargs="+",help="The pdbs")
    parser.add_argument('-n','--number',type=int,help='The number of pdbs to select')
    args = parser.parse_args()
    return args

def get_mp(args):
    molmap = {}
    for pdb in args.pdbs:
        print pdb
        pdbid = pdbtools.get_pdb_id(pdb)
        mpfilename = pdbid+"mpscore"
        mpfile = open(mpfilename,'r').readlines()
        for line in mpfile:
            if line.strip().startswith('MolProbity score'):
                score = line.split()[2]
                molmap[pdbid] = score
 
    return molmap

def select_top(args,molmap):
    if not os.path.isdir('topmp'):
        os.mkdir('topmp')
    count = 1
    for pdb in sorted(molmap.items(), key=operator.itemgetter(1),reverse=True):
        count+=1
        pdbid = pdb[0]
        for apdb in args.pdbs:
            if pdbtools.get_pdb_id(apdb) == pdbid:
                print pdbid
                os.system('cp %s topmp/'%apdb)
        if count > args.number:
            break
        
        

main()
