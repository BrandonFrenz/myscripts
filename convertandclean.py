#!/usr/bin/python
import argparse
import pdbtools
import amino_acids

def main():
    args = parseargs()
    clean_pdbs(args)
    

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdbs',nargs="+",help='The pdbs to clean')
    args = parser.parse_args()
    return args

def clean_pdbs(args):
    for pdb in args.pdbs:
        residues = pdbtools.get_unopened_residue_list(pdb)
        pdbtools.convert_to_standard_aas(residues)
        newresidues = []
        for res in residues:
            if res.name not in amino_acids.longer_names:
                continue
            else:
                newresidues.append(res)
        pdbtools.write_resis_to_pdb(newresidues,pdb)


main()
