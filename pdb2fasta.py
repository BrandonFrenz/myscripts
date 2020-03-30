#!/usr/bin/python
import argparse
import pdbtools
import amino_acids
import sys

def main():
    args = parseargs()
    print_sequence(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',help='The pdb')
    parser.add_argument('-sg','--structure_gaps',dest='structure_gaps',action='store_true',default=False)
    args = parser.parse_args()
    return args

def print_sequence(args):
    
    residues = pdbtools.get_unopened_residue_list(args.pdb)
    if args.structure_gaps:
        chains = []
        previousnum = 'x'
        
        sys.stdout.write('>'+pdbtools.get_pdb_id(args.pdb)+', '+str(len(residues)) +' residues\n')
        resit = 0
        sequence = []
        while resit < len(residues):
            residue = residues[resit]
            if args.structure_gaps and resit < len(residues) and resit != 0:
                if has_gap(residues[resit-1],residue) and sequence[-1] is not '/':
                    sys.stdout.write('/\n')
                    sequence.append('/')
            if (residue.chain not in chains and len(chains) != 0) or (previousnum != 'x' and residue.num != previousnum+1):
                if sequence[-1] != '/':
                    sys.stdout.write('/\n')
                    sequence.append('/')
            if residue.name in amino_acids.longer_names:
                sys.stdout.write(amino_acids.longer_names[residue.name])
                sequence.append(amino_acids.longer_names[residue.name])
            else:
                sys.stdout.write('X')
                sequence.append('/')
            previousnum = residue.num
            if residue.chain not in chains:
                chains.append(residue.chain)
            resit+=1
        sys.stdout.write('\n')
        #print sequence
    else:
        seqs = pdbtools.get_sequences(residues)
        for seq in seqs:
            print(''.join(seq))


def has_gap(residue1,residue2,dist_cutoff = 3.5):
    for atom in residue1.atoms:
        if atom.atomid == ' C  ':
            for atom2 in residue2.atoms:
                if atom2.atomid == ' N  ':
                    dist = pdbtools.atom_dist(atom,atom2)
                    if dist > dist_cutoff:
                        return True
                    else:
                        return False
    #this should not be reached
    print('No N and C atoms for the specificed residues were found')
    exit()

main()
