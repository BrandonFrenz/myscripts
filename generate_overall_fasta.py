#!/usr/bin/env python
import argparse
import pdbtools
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import amino_acids

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structures', nargs='+', help='The pdb files')
    parser.add_argument('-o', '--output', help='The output fasta')
    args = parser.parse_args()
    return args

def get_sequence(pdbfile):
    resis = pdbtools.get_unopened_residue_list(pdbfile)
    sequence = pdbtools.get_seq_from_resis(resis)
    return sequence

def combine_seq(seq1, seq2):
    alignment = pairwise2.align.localms(seq1, seq2, 2, 1, -.5, -.1)
    aln1 = alignment[0][0]
    aln2 = alignment[0][1]
    print(aln2)
    combined_seq = ''
    assert len(aln1) == len(aln2), 'The alignments arent the same length'
    for i in range(len(aln1)):
        aa1 = aln1[i]
        aa2 = aln2[i]
        if aa1 != '-' and aa2 != '-':
            print(aa1, aa2)
            assert aa1 == aa2, 'The sequences align but the amino acids dont match'
        if aa1 != '-':
            combined_seq+=aa1
        elif aa2 != '-':
            combined_seq+=aa2
    return combined_seq

def main():
    args = parseargs()
    sequence = get_sequence(args.structures[0])
    for i in range(1, len(args.structures)):
        pose_seq = get_sequence(args.structures[i])
        sequence = combine_seq(sequence, pose_seq)
    with open(args.output, 'w') as of:
        of.write('> combined sequence\n')
        of.write(sequence)



if __name__ == '__main__':
    main()
