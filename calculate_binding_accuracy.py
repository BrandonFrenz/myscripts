#!/usr/bin/python
import argparse
import os
import pdbtools
import amino_acids

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l','--ligand',help='The pdbfile of the ligand')
    parser.add_argument('-e','--experiment',nargs="+",help='The predicted ligands')
    parser.add_argument('-c', '--chains', nargs="+", default=[], help='Which chains are you interested in comparing')
    args = parser.parse_args()
    return args

def determine_accuracy( ligand, predictions, chains):
    ligands = [ligand]
    rms = determine_accuracies(ligands, predictions, chains)
    return rms

def determine_accuracies(ligands,predictions, chains):
    closest_rms = 1e9
    for ligand in ligands:
        ligand_resis = pdbtools.get_unopened_residue_list(ligand)
        for predicted in predictions:
            predicted_resis = pdbtools.get_unopened_residue_list(predicted)
            try:
                rms = get_rms(ligand_resis,predicted_resis, chains)
            except:
                print(ligand,predicted)
                print("EXITING")
                exit()
            if rms < closest_rms:
                closest_rms = rms
                print(predicted)
    return closest_rms


def get_rms(res1,res2,chains):
    backbones = ['N','C','CA','O', 'CB']
    atoms1 = []
    atoms2 = []
    for res in res1:
        if res.chain not in chains and len(chains) != 0:
            continue
        for atom in res.atoms:
            if atom.atomid.strip() in backbones or res.name not in amino_acids.longer_names:
                atoms1.append(atom)
    for res in res2:
        if res.chain not in chains and len(chains) != 0:
            continue
        for atom in res.atoms:
            if atom.atomid.strip() in backbones or res.name not in amino_acids.longer_names:
                atoms2.append(atom)
    rms = pdbtools.atomlist_rms(atoms1,atoms2)
    return rms

def main():
    args = parseargs()
    print(determine_accuracy(args.ligand, args.experiment, args.chains))

if __name__ == '__main__':
    main()
