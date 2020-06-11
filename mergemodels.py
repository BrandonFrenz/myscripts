#!/usr/bin/env python
import argparse
import pdbtools

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s1', '--structure1', help='The first structure')
    parser.add_argument('-s2', '--structure2', help='The second structure')
    parser.add_argument('-i1', '--insert1', type=int, help='The position to insert residues as')
    parser.add_argument('-cs', '--copy_start', type=int, help='The start of the residues to copy')
    parser.add_argument('-ce', '--copy_end', type=int, help='The end of the residues to copy')
    parser.add_argument('-o', '--output', help='The name of the output file')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    resis1 = pdbtools.get_unopened_residue_list(args.structure1)
    resis2 = pdbtools.get_unopened_residue_list(args.structure2)
    combined_resis = resis1[0:args.insert1-1]+resis2[args.copy_start-1:args.copy_end]+resis1[args.insert1-1:]
    pdbtools.write_resis_to_pdb(combined_resis, args.output, sort=False)

if __name__ == '__main__':
    main()
