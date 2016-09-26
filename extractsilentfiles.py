#!/usr/bin/python
import glob
import os
import multiprocessing
import sys
import argparse
import re

def main():
    args = parseargs()
    if args.mode == 'all':
        extractsilentfiles(args)
    if args.mode == 'percentage' or args.mode == 'topN':
        allmodels = readsilentfiles(args)
        tagmap = mapsilentfile(args,allmodels)
        selectedtags = toptags(args,tagmap)
        grabandwrite(allmodels,selectedtags)
        #if args.extract == True:
        extractfile(args,'selected.silent',1)
    if args.mode == 'count':
        allmodels = readsilentfiles(args)
    if args.mode == 'si':
        print_score_index(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m','--mode',help='which pdbs do you want to extract. All, top percentage, or top N, use si to print the score indexing')
    parser.add_argument('-s','--silentfiles',nargs="+",help='the silent files to extract')
    parser.add_argument('-n','--number',type=float,help='the percentage cutoff or number of poses to extract')
    parser.add_argument('-c','--cores',type=int,default=1,help='the number of computers to use')
    parser.add_argument('-si','--scoreindex',type=int,default=1,help='The, 0 base, indexing for the score you want to use in the silentfile')
    parser.add_argument('-ex','--executable',default='~/Desktop/Rosetta/main/source/bin/extract_pdbs.default.linuxgccrelease',help='The path to the extract pdbs executable')
    args = parser.parse_args()
    return args

#Change this command to point to your executable of the extract_pdbs app
def extractfile(args,silentfile,jobcount):
    print silentfile+" "+str(jobcount)
    path = "~/Desktop/Rosetta/main/source/bin/extract_pdbs.default.linuxgccrelease"
    os.system('%s -in:file:silent %s -out::prefix %d -crystal_refine -inout:skip_connect_info' % (args.executable,silentfile,jobcount))

def toptags(args,tagmap):
    totalposes = len(sorted(tagmap.items(), key=lambda x:x[1],reverse=True))
    tagcutoff = 0
    print totalposes
    if args.mode == 'topN':
        #all values less than tagcutoff or ignored
        tagcutoff = totalposes-args.number
    if args.mode == 'percentage':
        tagcutoff = totalposes*(1-(args.number/100))
    targettags = []
    count = 0
    for i in sorted(tagmap.items(), key=lambda x:x[1],reverse=True):
        if count < tagcutoff:
            count+=1
            continue
        targettags.append(i[0])
        count+=1
    return targettags


def readsilentfiles(args):
    allstructures = []
    header = True
    sfilecount = 1
    totalpdbs = 0
    for sfile in args.silentfiles:
        with open(sfile,'r') as silentfile:
            linecount = 1
            tag = 'thisisatemptag'
            newtag = 'thisisalsoatemptag'
            tagcounter = 1
            for line in silentfile:
                if linecount <= 3 and header == True:
                    allstructures.append(line)
                elif linecount > 3:
                    if line.startswith('SCORE:'):
                        totalpdbs +=1
                        tag = line.split()[-1]
                        newtag = ('tag_%s_%s'%(sfilecount,tagcounter))
                        tagcounter+=1
                    line = re.sub(tag,newtag,line)
                    allstructures.append(line)
                    header = False
                linecount+=1
            sfilecount +=1
    print " read " + str(totalpdbs) + " pdbs from the silentfiles"
    return allstructures

def mapsilentfile(args,silentfile):
    tagmap = {}
    firstscore = True
    for line in silentfile:
        if line.startswith('SCORE:'):
            if firstscore == False:
                data = line.split()
                score = float(data[args.scoreindex])
                tag = data[-1]
                tagmap[tag] = score
            if firstscore == True:
                data = line.split()
                print 'sorting by the lowest ',data[args.scoreindex],'value'
            firstscore = False
    return tagmap

#grab the structures matching the provided tags from the provided silentfile and write the silentfile to disk
def grabandwrite(silentfile,tags):
    selectedmodels = []
    header = True
    linecount = 1
    addstruct = False
    for line in silentfile:
        if header == True:
            selectedmodels.append(line)
        else:
            if line.startswith('SCORE:'):
                data = line.split()
                tag = data[-1]
                if tag in tags:
                    addstruct = True
                else:
                    addstruct = False
            if addstruct == True:
                selectedmodels.append(line)
        linecount+=1
        if linecount > 3:
            header = False
    with open('selected.silent','w') as silentfile:
        for line in selectedmodels:
            silentfile.write(line)

def print_score_index(args):
    with open(args.silentfiles[0],'r') as sf:
        for line in sf:
            if line.startswith('SCORE:'):
                data = line.split()
                it = 0
                for dat in data:
                    if it != 0:
                        print it,dat
                    it+=1
                break;

def extractsilentfiles(args):
    args.cores
    jobcount = 0
    silentfiles = args.silentfiles
    totaljobs = len(args.silentfiles)
    while jobcount < totaljobs:
        if __name__ == '__main__':
            jobs = []
            for i in range(0,args.cores):
                silentfile = silentfiles[jobcount]
                p = multiprocessing.Process(target=extractfile, args=(args,silentfile,jobcount))
                jobs.append(p)
                p.start() 
                jobcount = jobcount + 1
                if jobcount >= totaljobs:
                    break
            for i in jobs:
                i.join()

main()
