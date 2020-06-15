#!/usr/bin/env python
import argparse
import pdbtools

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure', help='The structure file')
    parser.add_argument('-o', '--output', help='The output file')
    parser.add_argument('-r', '--residues', type=int, nargs='+', help='The residues to convert')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    residues = pdbtools.get_unopened_residue_list(args.structure)
    newresis = []
    for res in residues:
        if res.num in args.residues:
            res = pdbtools.convert_to_gly(res)
        newresis.append(res)
    pdbtools.write_resis_to_pdb(newresis, args.output, False)

if __name__ == '__main__':
    main()
