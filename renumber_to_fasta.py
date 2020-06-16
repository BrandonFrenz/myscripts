#!/usr/bin/env python
import argparse
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import pdbtools

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure', help='The structure file')
    parser.add_argument('-c', '--chain', help='The relevant chain')
    parser.add_argument('-f', '--fasta', help='The fasta file')
    parser.add_argument('-o', '--output', help='The name of the output file')
    parser.add_argument('-os', '--offset', type=int, default=0, help='Use this to add missing numbers to the nterm, i.e. fasta starts at res 2')
    parser.add_argument('-kh', '--keep_het', defualt=False, action="store_true", help='Add all the het atoms back to the residue list at the end')
    args = parser.parse_args()
    return args

def align(pdbseq, fastaseq):
    alignments = pairwise2.align.localms(pdbseq, fastaseq, 2, 1, -.5, -.1)
    return alignments[0][0]

def parse_fasta(fastafile):
    seq = ''
    with open(fastafile,'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                continue
            line = line.strip()
            seq+=line
    return seq

def renumber_residues(resis, chain, aln_seq, offset):
    seqnum = offset
    resindex = 0
    for aa in aln_seq:
        seqnum+=1
        if aa != '-':
            resis[resindex].num = seqnum
            resindex+=1
    return resis

def main():
    args = parseargs()
    fasta_seq = parse_fasta(args.fasta)

    resis = pdbtools.get_unopened_residue_list(args.structure)
    target_resis = pdbtools.get_chain_resis(resis, args.chain)
    seq = ''.join(pdbtools.get_sequence(target_resis, False))
    
    aln_seq = align(seq, fasta_seq)

    renumbered_resis = renumber_residues(target_resis, args.chain, aln_seq, args.offset)
    if args.keep_het == True:
        het_resis = []
        for res in resis:
            record = res.atoms[0].record
            if record == 'HETATM':
                het_resis.append(res)
        for res in het_resis:
            renumbered_resis.append(res)
    pdbtools.write_resis_to_pdb(renumbered_resis, args.output, False)

if __name__ == '__main__':
    main()
