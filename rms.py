#!/usr/bin/python
import argparse
import pdbtools


#calculates the RMSD over the backbone atoms. Assumes structures are equal length
def main():
    args = parseargs()
    print_rms(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p1','--pdb1',help='the first pdb')
    parser.add_argument('-p2','--pdb2',help='the second pdb')
    args = parser.parse_args()
    return args

def print_rms(args):
    res1 = pdbtools.get_unopened_residue_list(args.pdb1)
    res2 = pdbtools.get_unopened_residue_list(args.pdb2)
    backbones = ['N','C','CA','O','CB']
    atoms1 = []
    atoms2 = []
    for res in res1:
        for atom in res.atoms:
            if atom.atomid.strip() in backbones:
                atoms1.append(atom)
    for res in res2:
        for atom in res.atoms:
            if atom.atomid.strip() in backbones:
                atoms2.append(atom)
    rms = pdbtools.atomlist_rms(atoms1,atoms2)
    print rms

main()
