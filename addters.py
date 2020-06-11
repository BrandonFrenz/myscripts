#!/usr/bin/env python
import argparse
import pdbtools

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structures', nargs='+', help='The structure files')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    for structure in args.structures:
        resis = pdbtools.get_unopened_residue_list(structure)
        for i in range(0, len(resis)-1):
            for atom in resis[i].atoms:
                if atom.atomid == ' C  ':
                    Catom = atom
                    break
            for atom in resis[i+1].atoms:
                if atom.atomid == ' N  ':
                    Natom = atom
                    break
            dist = pdbtools.atom_dist(Catom, Natom)
            if dist > 2:
                resis[i].isterm = True
        pdbtools.write_resis_to_pdb(resis, structure)


if __name__ == '__main__':
    main()
