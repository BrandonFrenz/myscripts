#!/usr/bin/python
import argparse
import pdbtools
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import amino_acids
import numpy

def main():
    args = parseargs()
    setup_cm(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--fasta',help='The fasta file')
    parser.add_argument('-p','--pdbs',nargs="+",help='The pdbs')
    args = parser.parse_args()
    return args

def setup_cm(args):
    alignment,mapped = find_unassigned_sequence(args)
    make_fasta(args,alignment,mapped)
    print_xml_lines(args)

#searches for any elements of the fasta file that are unmapped in the proteins
def find_unassigned_sequence(args):
    seq = ''
    with open(args.fasta,'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                continue
            line = line.strip()
            seq+=line
    unmapped = [False]*len(seq)
    for pdb in args.pdbs:
        resis = pdbtools.get_unopened_residue_list(pdb)
        pdbseq = pdbtools.get_seq_from_resis(resis)
        alignments = pairwise2.align.localms(pdbseq,seq,2,1,-.5,-.1)
        alignstring = alignments[0][0]
        start = ''
        startcount = 0
        end = ''
        endcount = 0
        previousaa = ''
        done = False
        while not done:
            aacount = 0
            done = True
            for aa in alignstring:
                if aa != '-' and previousaa not in amino_acids.amino_acids:
                    start = aa
                    startcount = aacount
                    if start == end:
                        #uses the format XXXEND----STARTXXX
                        print start,end
                        assert(end < start), 'the end of the first gap is not less than the start of the second'
                        alignstring = resolve_conflict(args,alignstring,res1,res2,endcount,startcount)
                        done = False
                elif previousaa != '-' and aa == '-':
                    end = previousaa
                    endcount = aacount
                previousaa = aa
                aacount+=1
        for i in range(0,len(alignstring)):
            if alignstring[i] in amino_acids.amino_acids:
                unmapped[i] = True
    return seq,unmapped

def resolve_conflicts(args,alignstring,res1,res2,endcount,startcount):
    res1 = res1[endcount]
    N = ''
    for atom in res1.atoms:
        if atom.name.strip() == 'N':
            N = atom
    res2 = res1[startcount]
    C = ''
    for atom in res2.atoms:
        if atom.name.strip() == 'C':
            C = atom
    dist = pdbtools.atom_dist(N,C) 

def make_fasta(args,seq,mapped):
    with open('hybrid.fasta','w') as fastafile:
        fastafile.write('>This is a custom fasta for hybrid with %s residues\n'%(len(seq)))
        newstr = ''
        first_seq = True
        for i in range(0,len(seq)):
            if mapped[i] == True:
                newstr+=seq[i]
                if i != len(seq):
                    if mapped[i] == False:
                        newstr+='/\n'
                firsts_seq = False
        fastafile.write(newstr)

def print_xml_lines(args):
    for pdb in args.pdbs:
        print '<Template pdb="%s" weight="1.0" cst_file="AUTO"/>'%pdb

main()
