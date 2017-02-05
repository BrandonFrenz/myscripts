#!/usr/bin/python
import argparse
from operator import itemgetter
import os

def main():
    args = parseargs()
    get_pdbs(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-hp','--hyak_path',help='The path on hyak')
    parser.add_argument('-sc','--scorefile',help='The scorefile to use')
    parser.add_argument('-n','--number',type=int,help='The number of structures to download')
    args = parser.parse_args()
    return args

def get_pdbs(args):
    scorepairs = parse_scorefile(args)
    scorepairs = sorted(scorepairs, key=itemgetter(0))
    for i in range(0,args.number):
        pdbid = scorepairs[i][0]
        command = 'scp hyak:%s/%s.pdb ./'%(args.hyak_path,scorepairs[i][1])
        os.system(command)

def parse_scorefile(args):
    scorepairs = []
    with open(args.scorefile,'r') as sf:
        for line in sf:
            if line.startswith('SCORE:'):
                try:
                    data = line.split()
                    pdbid = data[-1]
                    score = float(data[1])
                    scorepairs.append((score,pdbid))
                except:
                    do_nothing = True
    return scorepairs

main()
