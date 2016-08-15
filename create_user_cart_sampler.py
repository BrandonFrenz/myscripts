#!/usr/bin/python
import argparse
import amino_acids
import fileinput
import sys
import re

def main():
    args = parseargs()
    badres = get_residues(args)
    print badres
    add_bad_residues_to_xml(args,badres)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-x','--templatexml',help='The template xml')
    parser.add_argument('-r','--residues',help='A test file with the list of residues to use')
    parser.add_argument('-a','--add',type=int,default=0,help='The amount to change the input residues by (use this to adjust to rosetta numbering')
    args = parser.parse_args()
    return args

def get_residues(args):
    resfile = open(args.residues,'r').readlines()
    residues = []
    for line in resfile:
        line = line.split()
        lower = int(line[0])
        upper = lower
        if len(line) > 1:
            upper = int(line[1])
        for i in range(lower,upper+1):
            residues.append(str(i+args.add))
    return residues


def add_bad_residues_to_xml(args,badresidues):
    resline = 'residues='+','.join(badresidues)
    for line in fileinput.input(args.templatexml,inplace=1):
        if 'strategy="user"' in line:
            line = re.sub('residues=\d+',resline,line)
        sys.stdout.write(line)


main()
