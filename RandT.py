#!/usr/bin/python
#this script takes a protein dimer in pdb format and performs a random rotation and translation on one chain
import argparse
import pdbtools
import numpy as np
import math
import random

def main():
    args = parseargs()
    apply(args)


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',help='The pdbfile to modify')
    parser.add_argument('-o','--output',default='R&T.pdb',help='The name of the output file')
    parser.add_argument('-c','--chains',default=['B'],nargs='+',help='The chains to move')
    args = parser.parse_args()
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
    #axis = pdbtools.get_center_of_mass(resis)
    axis = [random.uniform(1,359),random.uniform(1,359),random.uniform(1,359)]
    theta = random.uniform(1,359)
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis,axis))
    rot_matrix = get_rotation_matrix(axis,theta)
    
    #apply transform
    rot_resis = []
    for resi in resis:
        if resi.chain in args.chains:
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

def com_to_origin(args,residues):
    
    #get moving residues COM
    moving_resis = []
    for residue in residues:
        if residue.chain in args.chains:
            moving_resis.append(residue)
    com = pdbtools.get_center_of_mass(moving_resis)
    
    moved_residues = []
    for residue in residues:
        if residue.chain in args.chains:
            movedatoms = []
            for atom in residue.atoms:
                atom.x -= com[0]
                atom.y -= com[1]
                atom.z -= com[2]
                movedatoms.append(atom)
            residue.atoms = movedatoms
        moved_residues.append(residue)
    return moved_residues

def write_pdb(args,resis):
    pdblines = pdbtools.make_pdblines_from_residues(resis)
    with open(args.output,'w') as outfile:
        for line in pdblines:
            outfile.write(line)

main()
