#!/usr/bin/env python
import argparse
import shutil

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structures', nargs='+', help='The pdbfiles')
    parser.add_argument('-r', '--residues', help='The residues in question in format x-y')
    args = parser.parse_args()
    return args

def parse_residues(residues):
    start = int(residues.split('-')[0])
    stop = int(residues.split('-')[1])
    residues = []
    for i in range(start, stop+1):
        residues.append(i)
    return residues

def parse_res_scores(structure):
    inscores = False
    res_scores = {}
    with open(structure, 'r') as inf:
        for line in inf:
            if line.startswith('pose'):
                inscores = True
                continue
            if line.startswith('#END_POSE_ENERGIES'):
                inscores = False
            if inscores:
                data = line.split()
                score = float(data[-1])
                resnum = int(data[0].split('_')[-1])
                res_scores[resnum] = score
    return res_scores

def main():
    args = parseargs()
    residues = parse_residues(args.residues)
    res_ranks = []
    for structure in args.structures:
        res_scores = parse_res_scores(structure)
        range_score = 0
        for res in residues:
            if res in res_scores.keys():
                range_score+=res_scores[res]
        res_ranks.append((range_score, structure))
    for count, pair in enumerate(sorted(res_ranks)):
        score = pair[0]
        pdb = pair[1]
        rank=count+1
        newname = f'rank_{rank}_{score}.pdb'
        shutil.copy(pdb, newname)
        



if __name__ == '__main__':
    main()
