#!/usr/bin/env python
import argparse
import pdbtools
import amino_acids

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structures', nargs="+", help='The structures to clean')
    parser.add_argument('-l', '--het', default=False, action='store_true', help='Keep het atms and all other non protein atoms')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    for structure in args.structures:
        residues = pdbtools.get_unopened_residue_list(structure)
        if not args.het:
            clean_resis = []
            for res in residues:
                if res.name not in amino_acids.longer_names:
                    continue
                else:
                    clean_resis.append(res)
            residues = clean_resis
        pdbtools.write_resis_to_pdb(residues, structure)


if __name__ == '__main__':
    main()
