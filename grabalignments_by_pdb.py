#!/usr/bin/python
import argparse
import re
import pdbtools

def main():
    args = parseargs()
    make_alignment(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--alignment',help='The alignment file')
    parser.add_argument('-p','--pdbs',nargs="+",help='The pdbs you want')
    parser.add_argument('-o','--outfile',default='newalignment.aln',help='The name of the output')
    args = parser.parse_args()
    return args

def make_alignment(args):
    alignments = get_alignment_lines(args)
    pdbids = get_pdbids(args)
    with open(args.outfile,'w') as outfile:
        for pdbid in alignments.keys():
            if pdbid in pdbids:
                outfile.write(''.join(alignments[pdbid]))

def get_alignment_lines(args):
    alignments = {}
    with open(args.alignment,'r') as afile:
        alignment_lines = []
        first = True
        for line in afile:
            if line.startswith('##'):
                if not first:
                    alignments[pdbid] = alignment_lines
                    alignment_lines = []
                pdbid = parse_pdbid(line)
            alignment_lines.append(line)
            first = False
    return alignments

def parse_pdbid(line):
    pdbelement = line.split()[-1]
    pdbid = re.split('_',pdbelement)[0]
    return pdbid

def get_pdbids(args):
    pdbids = []
    for pdb in args.pdbs:
        pdbids.append(pdbtools.get_pdb_id(pdb))
    return pdbids



main()
