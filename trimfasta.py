#!/usr/bin/python
import argparse
import re

def main():
    args = parseargs()
    header,sequence = parse_fasta(args)
    if args.mode == 'count':
        print 'there are',len(sequence),'residues'
        exit()

    args = parse_ranges(args)
    sequence = trim_fasta(args,sequence)
    write_sequence(args,header,sequence)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--fasta',help='The fasta file')
    parser.add_argument('-o','--output',help='The output file')
    parser.add_argument('-r','--ranges',nargs="+",help='The ranges to keep. use formation x-y')
    parser.add_argument('-m','--mode',help='Use -m count to count the fasta')
    args = parser.parse_args()
    return args

def parse_ranges(args):
    numbers = []
    for numrange in args.ranges:
        nums = re.split('-',numrange)
        lower = int(nums[0])
        upper = int(nums[1])
        for i in range(lower,upper+1):
            numbers.append(i)
    args.numbers = numbers
    return args

def parse_fasta(args):
    header = ''
    sequence = []
    with open(args.fasta,'r') as fastafile:
        for line in fastafile:
            if line.startswith('>'):
                header = line
            else:
                line = list(line.strip())
                for aa in line:
                    sequence.append(aa)
    return header,sequence

def trim_fasta(args,sequence):
    newseq = []
    it = 0
    while it < len(sequence):
        if it+1 in args.numbers:
            newseq.append(sequence[it])
        it+=1
    return newseq

def write_sequence(args,header,sequence):
    with open(args.output,'w') as outfile:
        outfile.write(header)
        for ele in sequence:
            outfile.write(ele)

main()
