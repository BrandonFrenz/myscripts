#!/usr/local/bin/python2.7
import sys
import pdbtools
import amino_acids

def clean_pdbs(pdbs):
    for pdb in pdbs:
        residues = pdbtools.get_unopened_residue_list(pdb)
        clean_resis = []
        for residue in residues:
            if residue.name not in amino_acids.longer_names:
                continue
            else:
                clean_resis.append(residue)
        pdbtools.write_resis_to_pdb(clean_resis,pdb)
