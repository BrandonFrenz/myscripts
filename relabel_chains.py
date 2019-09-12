#!/usr/bin/env python
import argparse
import pdbtools

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure', help='The structure file. Pdb only')
    parser.add_argument('-o', '--output', help='The output structure')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    resis = pdbtools.get_unopened_residue_list(args.structure)
    relabeled = pdbtools.relabel_chains(resis)
    pdbtools.write_resis_to_pdb(relabeled, args.output)

if __name__ == '__main__':
    main()
