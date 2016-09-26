#!/usr/bin/python
import argparse
import pdbtools
import re

def main():
    args = parseargs()
    args = parse_frag_ranges(args)
    fragments = pdbtools.parse_fragments(open(args.fragfile,'r').readlines())
    if args.mode == 'filter':
        fragments = filter_frags(args,fragments)
        pdbtools.write_fragments(fragments,args.outfile)
    if args.mode == 'count':
        print 'there are',len(fragments),'residues that start fragments'

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-ff','--fragfile',help='The fragmentfile')
    parser.add_argument('-fr','--fragrange',nargs="+",help='The fragment ranges to keep. use formation x-y')
    parser.add_argument('-o','--outfile',default='filterfrags.txt',help='The output file')
    parser.add_argument('-m','--mode',default='filter',help='The mode of the script')
    args = parser.parse_args()
    return args

def parse_frag_ranges(args):
    fragranges = []
    for fragrange in args.fragrange:
        fragnums = re.split('-',fragrange)
        lower = int(fragnums[0])
        upper = int(fragnums[1])
        for i in range(lower,upper+1):
            fragranges.append(i)
    args.fragnums = fragranges
    return args

def filter_frags(args,fragments):
    filterFrags = []
    it=1
    for frag in fragments:
        if frag.num in args.fragnums:
            frag.num = it
            it+=1
            filterFrags.append(frag)
    return filterFrags

main()
