#!/usr/bin/python
import argparse
import re
from rosetta import *

def main():
    args = parseargs()
    args = parse_residues(args)
    print_torsions(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdbs',nargs='+',help='The pdbs to check')
    parser.add_argument('-r','--residues',nargs='+',help='The residues to check')
    args = parser.parse_args()
    return args

def parse_residues(args):
    residues = []
    for res in args.residues:
        if len(res.split('-')) == 1:
            residues.append(int(res))
        else:
            rr = re.split('-',res)
            lower = int(rr[0])
            upper = int(rr[1])
            for i in range(lower,upper+1):
                residues.append(i)
    args.residues = residues
    return args

def print_torsions(args):
    rosetta.init(extra_options = '-mute all')
    for pdb in args.pdbs:
        pose = rosetta.pose_from_pdb(pdb)
        for res in args.residues:
            print pdb, res, pose.phi(res), pose.psi(res), pose.omega(res)

main()
