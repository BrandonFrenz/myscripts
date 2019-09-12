#!/usr/bin/env python
import argparse
from pyrosetta import *
from pathlib import Path
init('-ignore_unrecognized_res')

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structures', nargs='+', help='The structure files')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    scorefxn = get_fa_scorefxn()
    scores = []
    for structure in args.structures:
        stem = Path(structure).stem
        pose = pose_from_file(structure)
        score = (scorefxn)(pose)
        scores.append((stem, score))
    for score in scores:
        print(score)

if __name__ == '__main__':
    main()
