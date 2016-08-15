#!/usr/bin/python
import argparse
import fileinput
import re

def main():
    args = parseargs()
    newfragfile = renumberfrags(args)
    writefragfile(args, newfragfile)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fragments', help='the fragment file to renumber')
    parser.add_argument('-s', '--start', type=int, help='the number to start at')
    parser.add_argument('-n', '--newfrags', default='renumbered.frags', help='the name you want for the new fragfile')
    args = parser.parse_args()
    return args

def renumberfrags(args):
    newfragfile = []
    newfragcounter = args.start
    with open(args.fragments,'r') as fragments:
        for line in fragments:
            if line.startswith('position'):
                splitline = line.split()
                currentpos = splitline[1]
                #keep the end of the number at the 13th space after position:
                stringtosub = " "*(13-len(currentpos))+currentpos
                newstring = " "*(13-len(str(newfragcounter)))+str(newfragcounter)
                line = re.sub(stringtosub,newstring,line,1)
                newfragcounter+=1

            newfragfile.append(line)
    return newfragfile

def writefragfile(args,fragfile):
    with open(args.newfrags,'w') as newfragfile:
        for line in fragfile:
            newfragfile.write(line)

main()
