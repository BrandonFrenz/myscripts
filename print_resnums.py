#!/usr/bin/python
import argparse
import pdbtools

def main():
    args = parseargs()
    residues = pdbtools.get_unopened_residue_list(args.pdb)
    if args.rosnum != None:
        print residues[args.rosnum-1].num,residues[args.rosnum-1].chain,residues[args.rosnum-1].name
    if args.pdbnum != None:
        count = 1
        for res in residues:
            if res.num == args.pdbnum:
                print count,res.chain,res.name
            count+=1

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',help='The pdbfile')
    parser.add_argument('-rnum','--rosnum',type=int,help='Print the residue that corresponds to this Rosetta number')
    parser.add_argument('-pnum','--pdbnum',type=int,help='Print the rosetta number of this pdb residue')
    args = parser.parse_args()
    return args

main()
