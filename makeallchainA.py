#!/usr/bin/env python
import argparse
import pdbtools

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure', help='The input structure')
    parser.add_argument('-o', '--output', help='The output file name')
    parser.add_argument('-c', '--chain', default='A', help='The target chain')
    parser.add_argument('-kh', '--keep_het', default=False, action='store_true', help='Leave the het atoms on the chain they started with')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    residues = pdbtools.get_unopened_residue_list(args.structure)
    for residue in residues:
        record = residue.atoms[0].record
        if args.keep_het and record == 'HETATM':
            continue
        residue.chain = args.chain
    pdbtools.write_resis_to_pdb(residues, args.output, sort=False)

if __name__ == '__main__':
    main()
