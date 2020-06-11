#!/usr/bin/env python
import argparse
import pdbtools

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure', help='The input structure')
    parser.add_argument('-o', '--output', help='The output file name')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    residues = pdbtools.get_unopened_residue_list(args.structure)
    for residue in residues:
        residue.chain = 'A'
    pdbtools.write_resis_to_pdb(residues, args.output, sort=False)

if __name__ == '__main__':
    main()
