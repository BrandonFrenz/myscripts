#!/usr/bin/python
import argparse
import re

def main():
    args = parseargs()
    rename_nn(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdbs',nargs="+",help='The pdbfiles')
    args = parser.parse_args()
    return args

def rename_nn(args):
    for pdb in args.pdbs:
        num = re.split('after_filter_nn|_',pdb)[2]
        command = 'cp %s nn_%s.pdb'%(pdb,num)
        print command

main()
