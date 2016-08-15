#!/usr/local/bin/python2.7
import sys
import re
import multiprocessing
import fileinput
import os
import argparse

def main():
    args = parseargs()
    if args.mode == "count":
        printfastacount(args)
    if args.mode == 'prepdock' or args.mode == 'rundock':
        fragdocking(args)
    if args.mode== 'prepassembly' or args.mode == 'runassembly':
        scoreandassemble(args)
    if args.mode == 'preppdb':
        prepmodel(args)
    if args.mode == 'prepgrower':
        prepgrower(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m','--mode',default='None', help='Tell the script what you want it to do? Options include, rundock, preppdb, etc...')
    parser.add_argument('-p','--pdb',default='round1_model.pdb', help='The pdb you want the script to operate on')
    parser.add_argument('-n','--name',default='input.pdb', help='The name of the pdb being produced.')
    parser.add_argument('-fc','--fragcutoff',default=12,type=int,help='All segments equal to or shorter than this value will be removed from the inputpdb')
    parser.add_argument("-r", "--residues", type=int, nargs="+", help="the residues to build")
    parser.add_argument("-c", "--cores", type=int, default=16, help="the number of cores you have available to run")
    parser.add_argument("-f", "--fasta", help="the fasta file for counting (regular running reads the file in 1A.sh)")
    args = parser.parse_args()
    return args

def printfastacount(args):
        rescount = count_residues(args.fasta)
        ntosearch = rescount*20
        n_filtered = rescount*10
        print str(rescount) + " total residues. Suggested params n_to_search " + str(ntosearch) + " n_filtered " + str(n_filtered)

def fragdocking(args):
    if args.mode == "rundock":
        runparalleldock(args)
    if argsmode == "prepdock":
        setupque(args)


def run(residue):
    os.system('sh 1A.sh %s' % residue)

def runparalleldock(args):
    residues = args.residues
    compnumber = args.cores
    jobcount = 0
    while jobcount < len(residues):
        if __name__ == '__main__':
            jobs = []
            it = 0
            while it < compnumber:
               if not os.path.isfile('round1/input.%s.silent'%residues[jobcount]):
                    p = multiprocessing.Process(target=run, args=(residues[jobcount],))
                    jobs.append(p)
                    p.start()
                    it+=1
               jobcount = jobcount + 1
               if jobcount >= len(residues):
                   break
            for x in jobs:
               x.join()

def setupque(args):
    args.cores
    fasta = ''
    with open('1A.sh', 'r') as runscript:
        for line in runscript:
            if 'in::file::fasta' in line:
                line = line.strip().split()
                fasta = line[1]
    os.system('mkdir round1')
    totalresidues = count_residues(fasta)
    fragsize = checkfragsize()
    totaljobs = totalresidues-fragsize
    newjobs = find_missing(totaljobs)
    totaljobs = len(newjobs)
    totalnodes = int(totaljobs / args.cores) + (totaljobs % args.cores > 0)
    lower = 0
    str1 = os.getcwd()
    str2=str1.split('/')
    n=len(str2)
    name = str2[n-1]
    for i in range(1, totalnodes+1):
        jobresidues = []
        upper = lower+args.cores
        if upper>totaljobs:
            upper = totaljobs
        for res in range(lower, upper):
            jobresidues.append(str(newjobs[res]))
        if len(jobresidues) == 0:
            break
        update_hyakrun(jobresidues,name, args.cores)
        os.system('qsub -q bf hyakrun.sh')
        lower = upper

def count_residues(fasta):
    with open(fasta, 'r') as fastafile:
        firstline = True
        totalresidues = 0
        for line in fastafile:
            if firstline == True:
                firstline = False
            else:
                line = line.replace('\\','')
                line = line.strip()
                totalresidues+=len(line)
        return totalresidues

def checkfragsize():
    with open('1A.sh','r') as runscript:
        for line in runscript:
            if '-fragfile' in line:
                line = line.strip().split()
                fragfile = line[1]
                with open(fragfile) as fragments:
                    robettastyle = False
                    firstspace = True
                    fragsize = 0
                    for line in fragments:
                        if line.startswith('FRAME'):
                            line = line.strip().split()
                            fragsize = int(line[2])-int(line[1])+1
                            return fragsize
                        if line.startswith('position'):
                            robettastyle = True
                            continue
                        if robettastyle == True and line in ['\n','\r\n'] and firstspace == True:
                            firstspace = False
                            continue
                        if robettastyle == True and line in ['\n','\r\n'] and firstspace == False:
                            return fragsize
                        if robettastyle == True and firstspace == False:
                            fragsize +=1
                    return fragsize


def update_hyakrun(residues, name, cores):
    walltime = "3:59:00"
    residuesarg = ' '.join(residues)
    for line in fileinput.input('hyakrun.sh', inplace=1):
        line = re.sub("-N \S+", "-N %s" % (name), line)
        line = re.sub("walltime=\S+", "walltime=%s" % (walltime), line)
        line = re.sub("ppn=\d+,", "ppn=%s," % (cores), line)
        line = re.sub("feature=\d+", "feature=%s" % (cores), line)
        if 'denovodensity' in line:
            line = 'python denovodensity.py -m rundock -r %s -c %s &> /dev/null \n' % (residuesarg, cores)
        cwd = os.getcwd()
        line = re.sub("-d \S+", "-d %s/" % (cwd), line)
        sys.stdout.write(line)

def find_missing(expectedtotal):
    missing = []
    for i in range(1, expectedtotal+1):
        if os.path.isfile('round1/input.%d.silent' % i) == False:
            missing.append(i)
    print missing
    return missing



def get_sequence(fasta):
    sequence = ''
    with open(fasta,'r') as inputfasta:
        first = True
        for line in inputfasta:
            if '>' not in line:
                line = line.strip()
                for i in line:
                    sequence+=i
            if not first and '>' in line:
                sequence+='/\n'
            first = False
    print sequence
    return sequence

def prepmodel(args):
    newpdb = []
    terminalres = []
    chain = []
    previousres = 1;
    chainstart = 1;
    chainend  = 0;
    chainseq = ''
    with open(args.pdb,'r') as dockerpdb:
        for line in dockerpdb:
            if line.startswith('ATOM'):
               numbers = [int(s) for s in line.split() if s.isdigit()]
               residue = numbers[1]
               if residue - previousres > 1:
                  chainend = previousres
                  if chainend-chainstart > args.fragcutoff and len(chain)>0:
                      for i in chain:
                          newpdb.append(i)
                  chain = []        
                  chainstart = residue
               #if residue != previousresidue:
                   #chainseq += 
               previousres = residue
               chain.append(line)
    with open(args.name,'w') as outputpdb:
        for i in newpdb:
            outputpdb.write(i)
    outputpdb.close()
    countmissingsegments(args,newpdb)

def countmissingsegments(args,pdbfile):
    segmentcount = 0
    previousres = 1
    for line in pdbfile:
        if line.startswith('ATOM'):
            numbers = [int(s) for s in line.split() if s.isdigit()]
            residue = numbers[1]
            if residue != previousres+1 and residue != previousres:
                segmentcount+=1
            previousres = residue
    countfasta = count_residues(args.fasta)
    if previousres != countfasta:
        segmentcount+=1
    return segmentcount

def scoreandassemble(args):
    if args.mode == 'prepsetup':
        prepsetup(args)
        os.system('qsub -q bf hyakassemble.sh')
    if args.mode == 'runassembly':
        runassembly(args)

def prepsetup(args):
    update_hyakassemble(args)

def update_hyakassemble(args):
    str1 = os.getcwd()
    str2=str1.split('/')
    n=len(str2)
    name = str2[n-1]
    walltime = 3
    for line in fileinput.input('hyakassemble.sh', inplace=1):
        line = re.sub("-N \S+", "-N %sS&A" % (name), line)
        line = re.sub("walltime=\S+", "walltime=3:59:00", line)
        line = re.sub("ppn=\d+,", "ppn=%s," % (cores), line)
        line = re.sub("feature=\d+", "feature=%s" % (cores), line)
        cwd = os.getcwd()
        line = re.sub("-d \S+", "-d %s/" % (cwd), line)
        sys.stdout.write(line)

def runassembly(args):
    os.system('sh 1B.sh')
    runparallelassembly(args.cores)
    os.system('sh 1D.sh')

def runsingleassembly(count):
    os.system('sh 1C.sh %s' % count)
    #print count

def runparalleljobs(compnumber):
    jobcount = 1
    while jobcount <= compnumber:
          if __name__ == '__main__':
             jobs = []
             for i in range(1, compnumber+1):
                 p = multiprocessing.Process(target=runsingleassembly, args=(jobcount,))
                 jobs.append(p)
                 p.start() 
                 jobcount = jobcount + 1
                 if jobcount > compnumber:
                     break
             for x in jobs:
                 x.join()

#Create directories for each loop in the missing pdb
def prepgrower(args):
    pdb = readpdb(args)
    missingsegments = countmissingsegments(args,pdb)
    it = 1
    cwd = os.getcwd()
    os.chdir('../')
    while it <= missingsegments:
        if not os.path.isdir(str(it)):
            os.mkdir(str(it))
        os.system('cp %s/* %s/'%(cwd,it))
        it+=1
    os.mkdir('comparator')
    
def readpdb(args):
    pdbfile = []
    with open(args.pdb) as pdb:
        for line in pdb:
            pdbfile.append(line)
    return pdbfile

main()
