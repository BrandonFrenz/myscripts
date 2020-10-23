#!/usr/bin/python
import argparse
import pdbtools

def main():
    args = parseargs()
    resis = pdbtools.get_unopened_residue_list(args.pdb)
    print len(resis)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',help='The pdbfile')
    args = parser.parse_args()
    return args

main()
