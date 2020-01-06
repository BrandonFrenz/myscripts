#!/usr/bin/python
import argparse
import pdbtools
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import amino_acids
import numpy
import time

###THIS SCRIPT DOESN"T WORK AND IS INCOMPLETE####

def main():
    args = parseargs()
    get_alignment(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p1','--pdb1',help='The first pdb')
    parser.add_argument('-p2','--pdb2',help='The second pdb')
    args = parser.parse_args()
    return args

def get_alignment(args):
    res1 = pdbtools.get_residue_list(open(args.pdb1,'r').readlines())
    res2 = pdbtools.get_residue_list(open(args.pdb2,'r').readlines())
    seq1 = pdbtools.get_seq_from_resis(res1)
    seq2 = pdbtools.get_seq_from_resis(res2)
    
    s = time.time()
    alignments = pairwise2.align.localms(seq1,seq2,2,1,-.5,-.1)
    e = time.time()
    print e-s
    #print alignments,len(alignments[0][0]),len(alignments[0][1])
    alignstring = alignments[0][0]
    start = ''
    startcount = 0
    end = ''
    endcount = 0
    previousaa = ''
    print alignstring
    done = False
    while not done:
        aacount = 0
        for aa in alignstring:
            if aa != '-' and previousaa not in amino_acids.amino_acids:
                start = aa
                startcount = aacount
                if start == end:
                    #uses the format XXXEND----STARTXXX
                    assert(end < start), 'the end of the first gap is not less than the start of the second'
                    align_string = resolve_conflict(alignstring,res1,res2,endcount,startcount)
                    break
                print 'start',start
            elif previousaa != '-' and aa == '-':
                end = previousaa
                endcount = aacount
                print 'end',end
            previousaa = aa
            aacount+=1

def resolve_conflicts(alignstring,res1,res2,endcount,startcount):
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

if '__name__' == __main__:
    main()
