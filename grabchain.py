#!/usr/bin/python
import argparse

def main():
    args = parseargs()
    grabchains(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',help='the pdb to grab')
    parser.add_argument('-c','--chains',nargs='+',help='the chains you want to grab')
    args = parser.parse_args()
    return args

def grabchains(args):
    with open(args.pdb,'r') as pdb:
        relevantTER = True
        for line in pdb:
            if line.startswith('HETNAM'):
                chain = line[15]
                if chain in args.chains:
                    print(line.strip())
                    relevantTER = True
            elif line.startswith('ATOM') or line.startswith('LINK'):
                chain = line[21]
                if chain in args.chains:
                    print(line.strip())
                    relevantTER = True
            elif line.startswith('SSBOND'):
                chain = line[15]
                if chain in args.chains:
                    print(line.strip())
                    relevantTER = True
            else:
                if line.startswith('TER'):
                    if relevantTER == True:
                        print(line.strip())
                relevantTER = False

main()
