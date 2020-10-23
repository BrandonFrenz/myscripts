#!/usr/bin/env python
import argparse
import shutil

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structures', nargs="+", help='The structures to rename')
    parser.add_argument('-rr', '--resrange', nargs="+", help='The resrange using format 1-N$chainid')
    args = parser.parse_args()
    return args

def parse_res_range(resranges):
    res_pairs = []
    for rr in resranges:
        data = rr.split('-')
        start = int(data[0])
        end = int(data[1])
        for i in range(start, end+1):
            res_pairs.append(i)
    return res_pairs

def parse_scores(pdbfile):
    scores = {}
    with open(pdbfile, 'r') as inf:
        inscores = False
        for line in inf:
            if line.startswith('pose'):
                inscores = True
                continue
            if line.startswith('#END_POSE'):
                inscores = False
                break
            if inscores:
                resnum = int(line.split()[0].split('_')[-1])
                score = float(line.split()[-1])
                scores[resnum] = score
    return scores

def main():
    args = parseargs()
    residues = parse_res_range(args.resrange)
    score_struct_pairs = []
    for structure in args.structures:
        scores = parse_scores(structure)
        resscore = 0
        for res in residues:
            resscore+=scores[res]
        score_struct_pairs.append((resscore, structure))
    for i, scorepair in enumerate(sorted(score_struct_pairs)):
        score = scorepair[0]
        newname = f'rank{i+1}_{score}.pdb'
        shutil.copy(scorepair[1], newname)
        print(newname)


if __name__ == '__main__':
    main()
