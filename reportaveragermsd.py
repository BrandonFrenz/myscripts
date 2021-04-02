#!/usr/bin/env python
import argparse
import pdbtools
import statistics

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structures', nargs="+", help='The structures to get the average residue rmsd of')
    parser.add_argument('-o', '--output', help='The output file name')
    args = parser.parse_args()
    return args

def weighted_average(weight, res1, res2):
    newatoms = []
    res2bbs = pdbtools.get_backbones(res2, False)
    for i, atom in enumerate(pdbtools.get_backbones(res1, False)):
        atom2 = res2bbs[i]
        atom.x = (weight*atom.x+atom2.x)/(weight+1)
        atom.y = (weight*atom.y+atom2.y)/(weight+1)
        atom.z = (weight*atom.z+atom2.z)/(weight+1)
        newatoms.append(atom)
    res1.atoms = newatoms



def get_average_backbone(structures):
    starting_resis = pdbtools.get_unopened_residue_list(structures[0])
    for i in range(1, len(structures)):
            next_resis = pdbtools.get_unopened_residue_list(structures[i])
            for j, res in enumerate(next_resis):
                weighted_average(i, starting_resis[j], res)
    return starting_resis

def get_rmsd_variability(average_bb, structures):
    rmsds = {}
    average_bb = pdbtools.strip_non_protein(average_bb)
    for structure in structures:
        resis = pdbtools.get_unopened_residue_list(structure)
        resis = pdbtools.strip_non_protein(resis)
        for i, res in enumerate(resis):
            if i not in rmsds.keys():
                rmsds[i] = []
            rmsd = pdbtools.backbone_rmsd(average_bb[i], res, False)
            rmsds[i].append(rmsd)
    return rmsds

    

def main():
    args = parseargs()
    average_bb = get_average_backbone(args.structures)
    rmsds = get_rmsd_variability(average_bb, args.structures)
    with open(args.output, 'w') as of:
        of.write('residue number (rosetta numbering), average rmsd from center\n')
        for res in rmsds.keys():
            meanrmsd = statistics.mean(rmsds[res])
            of.write(f'{res}, {meanrmsd}\n')
    

if __name__ == '__main__':
    main()
