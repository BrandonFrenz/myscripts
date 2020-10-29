#!/usr/bin/env python
import argparse
import pdbtools

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    resis = pdbtools.get_unopened_residue_list(args.structure)
    for res in resis:
        for atom in res.atoms:
            print(f'"{atom.occupancy}"')

if __name__ == '__main__':
    main()
