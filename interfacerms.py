#!/usr/bin/python
import argparse
import os
import glob
import json
import pdbtools
import amino_acids

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structures', nargs="+", help='The structures to analyze')
    parser.add_argument('-n', '--native', help='The native structure to compare the rms to')
    parser.add_argument('-d', '--distance', type=float, default=8, help='Interface distance')
    parser.add_argument('-o', '--output', default=None, help='If you want to write the results to a json file set this to a file namke')
    args = parser.parse_args()
    return args

def parse_pdbids(pdbid_file):
    with open(pdbid_file,'r') as pf:
        data = json.load(pf)
    return data

def find_interface_resis(receptor_resis, ligand_resis, minimum_distance):
    receptor_interface = []
    ligand_interface = []
    for rres in receptor_resis:
        for ratom in rres.atoms:
            for lres in ligand_resis:
                if rres in receptor_interface and lres in ligand_interface:
                    continue
                for latom in lres.atoms:
                    if rres in receptor_interface and lres in ligand_interface:
                        break
                    dist = pdbtools.atom_dist(ratom, latom)
                    if dist <= minimum_distance:
                        if rres not in receptor_interface:
                            receptor_interface.append(rres)
                        if lres not in ligand_interface:
                            ligand_interface.append(lres)
    return receptor_interface, ligand_interface

def get_rms(res1, res2, interface_resis):
    backbones = ['N','C','CA','O','CB']
    atoms1 = []
    atoms2 = []
    for res in res1:
        if res.name not in amino_acids.longer_names:
            continue
        if interface_resis != None and not is_interface(res, interface_resis):
            continue
        for atom in res.atoms:
            if atom.atomid.strip() in backbones or res.name not in amino_acids.longer_names:
                atoms1.append(atom)
    for res in res2:
        if res.name not in amino_acids.longer_names:
            continue
        if interface_resis != None and not is_interface(res, interface_resis):
            continue
        for atom in res.atoms:
            if atom.atomid.strip() in backbones or res.name not in amino_acids.longer_names:
                atoms2.append(atom)
    rms = pdbtools.atomlist_rms(atoms1,atoms2)
    return rms

def is_interface(res, interface_resis):
    for ires in interface_resis:
        if ires.name == res.name and ires.num == res.num:
            return True
    return False

def determine_accuracies(native, predictions, ligand_chain, interface_resis):
    closest_rms = 1e9
    closest_pdb = None
    results = {}
    native_ligand_resis = pdbtools.get_unopened_residue_list(native)
    native_ligand_resis = pdbtools.get_chain_resis(native_ligand_resis, ligand_chain)
    for predicted in predictions:
        predicted_resis = pdbtools.get_unopened_residue_list(predicted)
        predicted_resis = pdbtools.get_chain_resis(predicted_resis, ligand_chain)
        rms = get_rms(native_ligand_resis, predicted_resis, interface_resis)
        print(predicted, rms)
        results[predicted] = rms
        try:
            rms = get_rms(native_ligand_resis, predicted_resis, interface_resis)
        except:
            print(ligand,predicted)
            print("EXITING")
            exit()
        if rms < closest_rms:
            closest_rms = rms
            closest_pdb = predicted
    return closest_rms, closest_pdb, results

def get_ligand_chain(resis):
    starting_chain = resis[0].chain
    for res in resis:
        if res.chain != starting_chain:
            return res.chain

def main():
    args = parseargs()
    #Find interface residues
    results = {}
    resis = pdbtools.get_unopened_residue_list(args.native)
    ligand_chain = get_ligand_chain(resis)
    receptor_resis = pdbtools.get_chain_resis(resis, 'A')
    ligand_resis = pdbtools.get_chain_resis(resis, 'B')
    receptor_interface, ligand_interface = find_interface_resis(receptor_resis, ligand_resis, args.distance)
    print(len(receptor_interface), len(ligand_interface))
    print('interface found')

    rms, rms_pdb, rms_results = determine_accuracies(args.native, args.structures, ligand_chain = ligand_chain, interface_resis = None)
    interface_rms, interface_rms_pdb, interface_rms_results = determine_accuracies(args.native, args.structures, ligand_chain = ligand_chain, interface_resis = ligand_interface)
    results['global_rms'] = rms_results
    results['interface_rms'] = interface_rms_results
    if args.output != None:
        with open(args.output, 'w') as of:
            json.dump(results, of, indent=4)
    print('RMSD:', rms, rms_pdb, interface_rms, interface_rms_pdb)

if __name__ == '__main__':
    main()
