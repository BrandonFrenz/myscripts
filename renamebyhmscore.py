#!/usr/bin/env python
import argparse
import shutil
import os

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--scores', help='The score file')
    args = parser.parse_args()
    return args

def parse_scores(scorefile):
    scores = []
    with open(scorefile, 'r') as inf:
        for line in inf:
            print(line)
            tag = line.split(':')[0]
            score = float(line.split()[-1])
            scores.append((score, tag))
    return scores


def main():
    args = parseargs()
    scores = parse_scores(args.scores)
    scores = sorted(scores)
    for rank, pair in enumerate(scores):
        score = pair[0]
        count=rank+1
        newname = f'rank{count}_{score}.pdb'
        pdb = f'{pair[1]}'
        if not os.path.isfile(pdb):
            print(f'Skiipping {pdb}. File is missing')
            continue
        shutil.copy(pdb, newname)

if __name__ == '__main__':
    main()
