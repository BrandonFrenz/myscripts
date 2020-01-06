#!/usr/bin/python
import argparse
import pdbtools
import amino_acids

def main():
    args = parseargs()
    move_to_end(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',help='The pdbfile')
    parser.add_argument('-o','--output',help='The output file name')
    args = parser.parse_args()
    return args

def move_to_end(args):
    residues = pdbtools.get_unopened_residue_list(args.pdb)
    chains = pdbtools.get_chains(residues)
    newresis = []
    for chain in chains:
        chain_resis = pdbtools.get_chain_resis(residues,chain)
        newchainresis=[]
        for res in chain_resis:
            if res.name in amino_acids.longer_names:
                newchainresis.append(res)
        for res in chain_resis:
            if res.name not in amino_acids.longer_names:
                newchainresis.append(res)
        newresis+=newchainresis
    pdbtools.write_resis_to_pdb(newresis,args.output,False)

main()
