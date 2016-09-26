#!/usr/bin/python
import argparse
import pdbtools

def main():
    args = parseargs()
    residues = combine_pdbs(args)
    pdbfile = pdbtools.make_pdblines_from_residues(residues)
    pdbtools.write_pdb(pdbfile,args.output)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p1','--pdb1',help='The first pdb')
    parser.add_argument('-p2','--pdb2',help='The second pdb')
    parser.add_argument('-o','--output',default='combined.pdb',help='The name of the output pdb')
    args = parser.parse_args()
    return args

def combine_pdbs(args):
    res1s = pdbtools.get_residue_list(open(args.pdb1,'r').readlines())
    res2s = pdbtools.get_residue_list(open(args.pdb2,'r').readlines())
    lastnum = res1s[-1].num+1
    renumed = []
    for res in res2s:
        res.num = lastnum
        lastnum+=1
        renumed.append(res)
    total_residues = res1s+renumed
    return total_residues

main()
