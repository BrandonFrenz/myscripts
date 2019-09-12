#!/usr/bin/python
import argparse
import os
import shutil

def main():
    args = parseargs()
    scores = {}
    for sf in args.scorefiles:
        parse_scorefile(sf,scores)
    select_models(scores,args.ntoselect)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--scorefiles',nargs="+",help='The score files')
    parser.add_argument('-n','--ntoselect',type=int,default=20,help='The number of models to select')
    args = parser.parse_args()
    return args

def parse_scorefile(scorefile,scores):
    with open(scorefile,'r') as sf:
        for line in sf:
            if line.startswith('SCORE:'):
                try:
                    score = float(line.split()[1])
                    pdbid = line.split()[-1]+'.pdb'
                    scores[pdbid] = [score]
                except:
                    continue

def select_models(scores,n):
    if not os.path.isdir('selected'):
        os.mkdir('selected')
    count = 0
    for key, value in sorted(scores.iteritems(), key=lambda (k,v): (v,k)):
        if count == n:
            break
        if not os.path.isfile(key):
            print 'missing',key
        shutil.copy(key,'selected/')
        count+=1

main()
