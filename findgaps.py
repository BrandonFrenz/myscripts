#!/usr/bin/env python
import argparse
import pdbtools

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure', help='The structure of interest')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    residues = pdbtools.get_unopened_residue_list(args.structure)
    residues = pdbtools.strip_non_protein(residues)
    for i in range(0, len(residues)-1):
        res1 = residues[i]
        res2 = residues[i+1]
        gap = pdbtools.connection_distance(res1, res2)
        if gap > 1.5:
            print(res1.name, res1.num, res2.name, res2.num, gap)

if __name__ == '__main__':
    main()
