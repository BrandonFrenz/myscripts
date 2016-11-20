#!/usr/bin/python
import argparse
import pdbtools

def main():
    args = parseargs()
    get_max_ca_dist(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',help='The pdbfile')
    args = parser.parse_args()
    return args

def get_max_ca_dist(args):
    resis = pdbtools.get_unopened_residue_list(args.pdb)
    it = 0
    maxdist = 0
    while it < len(resis)-1:
        #print resis[it].num,resis[it+1].num
        if resis[it].num != resis[it+1].num-1:
            it+=1
            continue
        distance = pdbtools.atom_dist(pdbtools.get_ca(resis[it]),pdbtools.get_ca(resis[it+1]))
        if maxdist > distance or maxdist == 0:
            maxdist = distance
        it+=1
    print maxdist

main()
