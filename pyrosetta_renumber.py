#!/usr/bin/env python
import argparse
from pyrosetta import *
init('-ignore_unrecognized_res')

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure', help='The structure file')
    parser.add_argument('-o', '--outfile', help='The output file')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    pose = pose_from_file(args.structure)
    pose.pdb_info().obsolete(True)
    pose.dump_pdb(args.outfile)

if __name__ == '__main__':
    main()
