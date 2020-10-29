1#!/usr/bin/env python
import argparse
from pyrosetta import *
init('-ignore_unrecognized_res -ignore_zero_occupancy false')

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure', help='The structure file')
    parser.add_argument('-o', '--output', help='The output file')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    pose = pose_from_file(args.structure)
    pose.dump_pdb(args.output)

if __name__ == '__main__':
    main()
