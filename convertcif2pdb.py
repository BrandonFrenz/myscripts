#!/usr/bin/python
import argparse
from Bio.PDB.MMCIFParser import MMCIFParser
import os
from Bio.PDB.PDBIO import PDBIO

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--structure_files',nargs="+",help='The cif files to convert')
    args = parser.parse_args()
    return args

def convert_files(structure_files):
    for sf in structure_files:
        prefix = os.path.splitext(sf)[0]
        print prefix
        parser = MMCIFParser()
        structure = parser.get_structure(prefix,sf)
        io = PDBIO()
        io.set_structure(structure)
        io.save(prefix+'.pdb')

def main():
    args = parseargs()
    convert_files(args.structure_files)

if __name__ == '__main__':
    main()
