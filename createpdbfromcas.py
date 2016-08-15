#!/usr/bin/python
import argparse
import re
import os

def main():
    args = parseargs()
    writecas(args)


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--calphas',help='the calpha file produced by dock_into_density app')
    args = parser.parse_args()
    return args

def writecas(args):
    with open(args.calphas,'r') as calphas:
        for line in calphas:
            line = line.split()
            atomid = line[1]
            x = format(float(line[2]), '.3f')
            y = format(float(line[3]), '.3f')
            z = format(float(line[4]), '.3f')
            pdbline = list(" "*80)
            pdbline[0:5] = "ATOM  "
            pdbline[6:10] = " "*(5-len(atomid))+atomid
            pdbline[12:15] = " MG "
            pdbline[17:20] = "UNK"
            pdbline[21] = "A"
            pdbline[22:25] = " "*(4-len(atomid))+atomid
            pdbline[30:37] = FloatToPDBString(x)
            pdbline[38:45] = FloatToPDBString(y)
            pdbline[46:54] = FloatToPDBString(z)
            print "".join(pdbline)

def FloatToPDBString(number):
    numstr = str(number)
    numlist = numstr.split('.')
    intpart = numlist[0]
    decpart = numlist[1]
    pdbstr = " "*(4-len(intpart))+intpart+"."+decpart+" "*(3-len(decpart))
    return pdbstr

main()
os.system('rm selectedpoints.txt')
