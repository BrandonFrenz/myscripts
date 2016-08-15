#!/usr/bin/python
import argparse
import amino_acids
import fileinput
import sys
import re

def main():
    args = parseargs()
    badres = get_outliers(args)
    add_bad_residues_to_xml(args,badres)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mp','--mpscore',help='The MolProbity score file')
    parser.add_argument('-x','--templatexml',help='The template xml')
    parser.add_argument('-c','--chains',nargs="+",default="0",help='Grab only residues on these chains, default all')
    parser.add_argument('-as','--asymsize',type=int,default=0,help='The size of the asymetric unit')
    args = parser.parse_args()
    return args

def get_outliers(args):
    mpfile = open(args.mpscore,'r').readlines()
    bad_residues = []
    ramalines = False
    for line in mpfile:
        if '----------Ramachandran angles----------' in line:
            ramalines = True
        elif ('SUMMARY:' in line):
            ramalines = False
        if ramalines == True:
            if len(line.split()) < 1:
                continue
            resname = ''.join(list(line)[10:13])
            if resname in amino_acids.longer_names:
                chain = list(line)[3]
                resnum = int(''.join(list(line)[4:8]))
                if chain in args.chains or args.chains == "0":
                    bad_residues.append((chain,resnum))
    return bad_residues

def add_bad_residues_to_xml(args,badres):
    badresidues = []
    for res in badres:
        if args.asymsize != 0:
            if int(res[1]) < args.asymsize:
                badresidues.append(str(res[1]))
        else:
            badresidues.append(str(res[1]))
    resline = 'residues='+','.join(badresidues)
    for line in fileinput.input(args.templatexml,inplace=1):
        if 'strategy="user"' in line:
            line = re.sub('residues=\d+',resline,line)
        sys.stdout.write(line)


main()
