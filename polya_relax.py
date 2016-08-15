#!/usr/bin/python
import pdbtools
import argparse
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio.PDB import StructureAlignment
import Bio.PDB
import re
import os

def main():
    args = parseargs()
    args = parse_resranges(args)
    salign = make_alignment(args.pdbs[0],args.staticpdb)
    convert_all_to_ala(args,salign)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-sp','--staticpdb',help='The pdb to align all others to')
    parser.add_argument('-p','--pdbs',nargs="+",help='The pdbs you intend to modify. Note: They all must be single chain and the same sequence')
    parser.add_argument('-l','--ligands',help='The ligands to add to the pdbs after')
    parser.add_argument('-cr','--convert_res',nargs="+",help='The residues to convert to ala. The format is X-X and the numbers correspond to the static pdb')
    args = parser.parse_args()
    return args

def parse_resranges(args):
    resranges = []
    for resrange in args.convert_res:
        rr = re.split('-',resrange)
        lower = int(rr[0])
        upper = int(rr[1])
        resrange = [lower,upper]
        resranges.append(resrange)
    args.resranges = resranges
    return args

#this function will get the sequences from the pdbs and make a clustalw alignment.
def make_alignment(pdb1,pdb2):
    seq1 = str(get_sequence(pdb1))
    seq2 = str(get_sequence(pdb2))
    s1 = pdbtools.get_bio_structure(pdb1)
    pdbid = pdbtools.get_pdb_id(pdb1)
    testres = s1[0]['A'][(' ',475,' ')]
    s2 = pdbtools.get_bio_structure(pdb2)
    pdbtools.make_alignment_input(seq1,seq2)#produces an alignment file of the 2 sequences named alignment.fasta
    cline = ClustalwCommandline('clustalw',infile='alignment.fasta')
    cline()
    falign = AlignIO.read("alignment.aln", "clustal")
    salign = StructureAlignment(falign,s1,s2)
    return salign

def convert_all_to_ala(args,salign):
    s1 = pdbtools.get_bio_structure(args.staticpdb)
    resis_to_change = []
    for resrange in args.resranges:
        lower = resrange[0]
        upper = resrange[1]
        slower = 0
        supper = 0
        for ( r1, r2 ) in salign.get_iterator():
            if r1 is None or r2 is None:
                continue
            resid = r2.get_full_id()
            r2num = resid[3][1]
            if r2num == lower:
                slower = r1.get_full_id()[3][1]
            if r2num == upper:
                supper = r1.get_full_id()[3][1]
        if slower == 0 or supper == 0:
            assert ( slower > 0 and supper > 0 ), "There is no aligned residue in either your start or stop position"
        for i in range(slower,supper+1):
            resis_to_change.append(i)

    for pdb in args.pdbs:
        pdbfile = open(pdb,'r').readlines()
        pdbfile = pdbtools.convert_resis_to_ala(pdbfile,resis_to_change)
        pdbid = pdbtools.get_pdb_id(pdb)
        name = pdbid+'.pdb'
        write_pdb_file(pdbfile,name)

def write_pdb_file(pdbfile,name):
    with open(name,'w') as newfile:
        for line in pdbfile:
            newfile.write(line)

def convert_struct_to_ala(structure,resis_to_change):
    for model in structure:
        for chain in model:
            for res in chain:
                resid = res.id
                if resid[1] in resis_to_change:
                    alares = pdbtools.convert_bio_res(res,'ALA')
                    chain.detach_child(resid)
                    chain.attach(alares)
    return structure

def get_sequence(pdb):
    s = pdbtools.get_bio_structure(pdb)
    ppbuilder = PPBuilder()
    model_nr = 1
    poly_peps = ppbuilder.build_peptides(s,model_nr)
    for polypep in poly_peps:
        return polypep.get_sequence()

main()
