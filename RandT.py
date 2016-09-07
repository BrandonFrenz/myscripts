#!/usr/bin/python
#this script takes a protein dimer in pdb format and performs a random rotation and translation on one chain
import argparse
import pdbtools
import numpy as np
import math
import random
import amino_acids

def main():
    args = parseargs()
    if args.chains == ['second']:
        args = get_second_chain(args)
    apply(args)


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',help='The pdbfile to modify')
    parser.add_argument('-o','--output',default='R&T.pdb',help='The name of the output file')
    parser.add_argument('-c','--chains',default=['all'],nargs='+',help='The chains to move')
    args = parser.parse_args()
    return args

def get_second_chain(args):
    resis = pdbtools.get_residue_list(open(args.pdb,'r').readlines())
    firstchain = resis[1].chain
    for residue in resis:
        if residue.chain != firstchain:
            args.chains = [residue.chain]
            return args

def apply(args):
    resis = pdbtools.get_residue_list(open(args.pdb,'r').readlines())
    resis = rotate_residues(args,resis)
    resis = com_to_origin(args,resis)
    write_pdb(args,resis)

def get_rotation_matrix(axis,theta):
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
            [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
            [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def rotate_residues(args,resis):
    #get matrix
    chains = get_chains(args,resis)
    rot_resis = []
    for chain in chains:
        #adds all the non moving residues to the list
        if chain not in args.chains and args.chains[0] != 'all':
            for res in resis:
                if res.chain == chain:
                    rot_resis.append(res)
        #moves all the relevant residues
        else:
            print 'rotating',chain
            axis = [random.uniform(1,359),random.uniform(1,359),random.uniform(1,359)]
            theta = random.uniform(1,359)
            axis = np.asarray(axis)
            axis = axis/math.sqrt(np.dot(axis,axis))
            rot_matrix = get_rotation_matrix(axis,theta)
            
            #apply transform
            for resi in resis:
                if resi.chain == chain:
                    rot_atoms = []
                    for atom in resi.atoms:
                        xyz = [atom.x,atom.y,atom.z]
                        newcoords = np.dot(rot_matrix, xyz)
                        atom.x = newcoords[0]
                        atom.y = newcoords[1]
                        atom.z = newcoords[2]
                        rot_atoms.append(atom)
                    resi.atoms = rot_atoms
                    rot_resis.append(resi)

    return rot_resis

#move each chain so that it's COM is at the origin
def com_to_origin(args,residues):
    
    chains = get_chains(args,residues)
    moved_residues = []
    for chain in chains:
        if chain not in args.chains and args.chains[0] != 'all':
            continue
        #get moving residues COM
        moving_resis = []
        for residue in residues:
            if residue.chain == chain:
                moving_resis.append(residue)
        com = pdbtools.get_center_of_mass(moving_resis)
        
        for residue in moving_resis:
            if residue.chain == chain:
                movedatoms = []
                for atom in residue.atoms:
                    atom.x -= com[0]
                    atom.y -= com[1]
                    atom.z -= com[2]
                    movedatoms.append(atom)
                residue.atoms = movedatoms
                moved_residues.append(residue)
    #add unmoved residues back to the list
    for chain in chains:
        if chain not in args.chains and args.chains[0] != 'all':
            for residue in residues:
                if residue.chain == chain:
                    moved_residues.append(residue)
    return moved_residues

def get_chains(args,resis):
    chains = []
    for resi in resis:
        if resi.name not in amino_acids.longer_names:
            continue
        if resi.chain not in chains:
            chains.append(resi.chain)
    return chains
    

def write_pdb(args,resis):
    pdblines = pdbtools.make_pdblines_from_residues(resis)
    with open(args.output,'w') as outfile:
        for line in pdblines:
            outfile.write(line)

main()
