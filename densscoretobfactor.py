#!/usr/bin/env python
import argparse
import pdbtools

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure', help='The pdb file of the structure')
    parser.add_argument('-o', '--output', help='The output file')
    args = parser.parse_args()
    return args

def parse_dens_scores(structure):
    scoresection = False
    dens_scores = []
    with open(structure, 'r') as inf:
        for line in inf:
            if line.startswith('label'):
                densindex = line.split().index('elec_dens_fast')
            if line.startswith('pose'):
                scoresection=True
                continue
            if line.startswith('#END_POSE'):
                scoresection=False
            if scoresection:
                if line.startswith("VRT"):
                    continue
                dens_score = line.split()[densindex]
                dens_scores.append(float(dens_score))
    return dens_scores


def main():
    args = parseargs()
    dens_scores = parse_dens_scores(args.structure)
    residues = pdbtools.get_unopened_residue_list(args.structure)
    for i, res in enumerate(residues):
        dens_score = dens_scores[i]
        for atom in res.atoms:
            temp = str(round(dens_score, 2))
            atom.tempfact = str(round(dens_score, 2))
    pdbtools.write_resis_to_pdb(residues, args.output)

if __name__ == '__main__':
    main()
