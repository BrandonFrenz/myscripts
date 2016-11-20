#!/usr/bin/python
import argparse
import re
import os

def main():
    args = parseargs()
    invert_rms_and_scores(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdbs',nargs="+",help='The pdbs')
    args = parser.parse_args()
    return args

def invert_rms_and_scores(args):
    for pdb in args.pdbs:
        'skip nn pdb'
        if 'nn' in pdb:
            continue
        print pdb
        cycle = re.split('_',pdb)[2]
        first = re.split('_',pdb)[4]
        second = re.split('_',pdb)[5]
        second = re.sub('.pdb','',second)
        newname = 'after_filter_'+cycle+'__'+second+'_'+first+'.pdb'
        os.system('mv %s %s'%(pdb,newname))

main()
