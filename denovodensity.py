#!/usr/local/bin/python2.7
import sys
import glob
import re
import math
import multiprocessing
import os
import fileinput
import argparse
import time

def main():
    args = parseargs()
    #The Setup
    if args.mode == 'prep':
        if args.step == 'full':
            runfull(args)
        elif args.step == 'dock':
            setupdocking(args)
        elif args.step == 'assemble':
            if not os.path.isfile('round%s/scores1'%(args.roundcount)):
                print " did not find scores1 in round%s/ folder generate it with step = score "%args.roundcount
                exit()
            queassembly(args)
        elif args.step == 'growerprep':
            grower_prep(args)
        #This block handles both scoring and consensus (They only need 1 core)
        elif args.step == 'pickfrags':
            que_fragpicking(args)
        else:
            prepsetup(args)
            os.system('qsub -q bf hyakrun.sh')
    #The running
    elif args.mode == 'run':
        if args.step == 'dock':
            runparalleldocking(args)
        elif args.step == 'score':
            os.system('sh 1B.sh')
            if args.track:
                updatetrackfile()
        elif args.step == 'assemble':
            runparallelassembly(args)
        elif args.step == 'consensus':
            os.system('sh 1D.sh')
            if args.track:
                updatetrackfile()
        elif args.step == 'pickfrags':
            run_fragpicking(args)


#Read the command line arguments
def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m','--mode',default='prep',help='The mode to use. Modes include prep and run')
    parser.add_argument('-n','--nodes',type=int,default=13,help='The number of nodes to use in step C. The default 13 assumes you are running nstruct 5')
    parser.add_argument('-c','--cores',type=int,default=16,help='The number of cores to use')
    parser.add_argument('-sc','--startcount',type=int,default=1,help='Where to start the counting for the output')
    parser.add_argument("-r", "--residues", type=int, nargs="+", help="the residues to dock")
    parser.add_argument('-s','--step',help='Which step of the protocol are you on?')
    parser.add_argument('-name','--name',help='what to name the job')
    parser.add_argument('-t','--track',type=bool,default=True,help='Track the jobs')
    parser.add_argument('-rc','--roundcount',type=int,default=1,help='Which round of docking you are on')
    parser.add_argument('-fs','--fragsize',type=int,default=20,help='The size cut off for fragment into the grower')
    parser.add_argument('-tr','--trim',type=int,default=3,help='The number of residues to trim off the starting inputs')
    parser.add_argument('-f','--fasta',default='full.fasta',help='The fasta file to use')
    args = parser.parse_args()
    return args

#update the hyak submission script for a single job
def prepsetup(args):
    update_hyakrun(args,args.startcount)

#update the hyakrun script with step B-D params
def update_hyakrun(args,count):
    str1 = os.getcwd()
    str2=str1.split('/')
    n=len(str2)
    for line in fileinput.input('hyakrun.sh', inplace=1):
        line = re.sub("-N \S+", "-N %s" % (args.name), line)
        line = re.sub("walltime=\S+", "walltime=3:59:00", line)
        line = re.sub("ppn=\d+,", "ppn=%s," % (args.cores), line)
        line = re.sub("feature=\d+", "feature=%s" % (args.cores), line)
        cwd = os.getcwd()
        line = re.sub("-d \S+", "-d %s/" % (cwd), line)
        if 'denovodensity.py' in line:
            line = 'python denovodensity.py -m run -sc %s -c %s -s %s -t %s -rc %s &> /dev/null \n'%(count,args.cores,args.step,args.track,args.roundcount)
        sys.stdout.write(line)

#update the hyakrun script with docking arguments
def update_hyakrun_dock(args,residues):
    walltime = "3:59:00"
    residuesarg = ' '.join(residues)
    for line in fileinput.input('hyakrun.sh', inplace=1):
        line = re.sub("-N \S+", "-N %s" % (args.name), line)
        line = re.sub("walltime=\S+", "walltime=%s" % (walltime), line)
        line = re.sub("ppn=\d+,", "ppn=%s," % (args.cores), line)
        line = re.sub("feature=\d+", "feature=%s" % (args.cores), line)
        if 'denovodensity' in line:
            line = 'python denovodensity.py -m run -s dock -r %s -c %s -t %s -rc %s &> /dev/null \n' % (residuesarg,args.cores,args.track,args.roundcount)
        cwd = os.getcwd()
        line = re.sub("-d \S+", "-d %s/" % (cwd), line)
        sys.stdout.write(line)

#runs the full score and assembly on a single node (SLOW).
def runfull(args):
    done = False
    args.step = steptrack()
    while not done:
        print "looping for another round"
        if args.step == 'dock':
            setupdocking(args)
            wait_for_jobs()
            args.step = 'score'
        if args.step == 'score':
            createtrackfile('score',1)
            prepsetup(args)
            os.system('qsub -q bf hyakrun.sh')
            wait_for_jobs()
            args.step = 'assemble'
        if args.step == 'assemble':
            queassembly(args)
            wait_for_jobs()
            args.step = 'consensus'
        if args.step == 'consensus':
            createtrackfile('consensus',1)
            prepsetup(args)
            os.system('qsub -q bf hyakrun.sh')
            wait_for_jobs()
            args.step = 'verify'
        if args.step == 'verify':
            percent = percent_finished(args)
            print str(percent) + " complete"
            if percent > 70:
                done = True
            newresidues = total_added(args)
            print str(newresidues) + " new residues added in the last round"
            if newresidues < 5:
                done = True
            if not done:
                args.step = 'dock'
                os.system('cp round%s_model.pdb round%sinput.pdb'%(args.roundcount,args.roundcount+1))
                args.roundcount+=1
                update_submissions(args)
    args.step = 'pickfrags'
    createtrackfile('pickfrags',1)
    update_hyakrun(args,0)
    os.system('qsub -q bf hyakrun.sh')
    wait_for_jobs()
    grower_prep(args)

#create the file for jobtracking
def createtrackfile(step,total):
    with open('jobtracking.txt','w') as trackfile:
        trackfile.write(str(step) + " 0/%s \n"%(total))

#update the trackfile by 1
def updatetrackfile():
    if not os.path.isfile('jobtracking.txt'):
        return
    for line in fileinput.input('jobtracking.txt', inplace=1):
        data = re.split(' |/',line)
        jobcount = int(data[1])
        jobcount +=1
        data[1]
        line = "%s %s/%s"%(data[0],jobcount,data[2])
        sys.stdout.write(line)

#check if all the jobs in the current step are finished
def checktrackfile():
    with open('jobtracking.txt','r') as trackfile:
        for line in trackfile:
            data = re.split(' |/',line)
            #You can stop jobs by changing the first word to be "done" in the trackfile
            if data[0] == 'done':
                exit()
            currentjob = data[1]
            totaljobs = data[2]
            if currentjob == totaljobs:
                return True
            else:
                return False

#return the step the current protocol is on according to the track file.
def steptrack():
    steps = ['dock','score','assemble','consensus','verify']
    step = 'dock'
    if os.path.isfile('jobtracking.txt'):
        with open('jobtracking.txt','r') as trackfile:
            for line in trackfile:
                split = line.split()
                step = split[0]
                split = re.split(' |/',line)
                current = split[1]
                total = split[2]
                #if the current step is not finished exit the protocol
                if int(current) != int(total):
                    print "Warning the current trackfile believes there are jobs still running"
                    exit()
            index = steps.index(step)
            #find the step after the current
            index+=1
            if index >= len(steps):
                index = 0
            step = steps[index]
    return step

#check the trackfile every minute until all jobs are finished
def wait_for_jobs():
    minutes = 0
    done = False
    while not done:
        done = checktrackfile()
        time.sleep(60)
        minutes+=1
        print " waiting for jobs " + str(minutes) + " minutes elapsed"
        if minutes > 6000:
            createtrackfile('timelimit',6000)
            exit()

def update_submissions(args):
    for line in fileinput.input('1A.sh', inplace=1):
        if "-out:file:silent" in line:
            line = "\t-out:file:silent round%s/input.$1.silent \\\n"%(args.roundcount)
        if args.roundcount > 1:
            if '-startmodel' in line:
                line = "\t-startmodel round%sinput.pdb \\\n"%args.roundcount
        sys.stdout.write(line)
    for line in fileinput.input('1B.sh', inplace=1):
        if "-in::file::silent" in line:
            line = "\t-in::file::silent round%s/input*silent \\\n"%args.roundcount
        if "-scorefile" in line:
            line = "\t-scorefile round%s/scores1 \\\n"%args.roundcount
        sys.stdout.write(line)
    for line in fileinput.input('1C.sh', inplace=1):
        if "-in::file::silent" in line:
            line = "\t-in::file::silent round%s/input*silent \\\n"%args.roundcount
        if "-scorefile" in line:
            line = "\t-scorefile round%s/scores1 \\\n"%args.roundcount
        if "-out:file:silent" in line:
            line = "\t-out:file:silent round%s/assembled.$1.silent \\\n"%(args.roundcount)
        sys.stdout.write(line)
    for line in fileinput.input('1D.sh', inplace=1):
        if "-in::file::silent" in line:
            line = "\t-in::file::silent round%s/assembled*silent \\\n"%args.roundcount
        if "cp S_0001.pdb" in line:
            line = "cp S_0001.pdb round%s_model.pdb \n"%args.roundcount
        sys.stdout.write(line)


#submit the assembly on multiple nodes
def queassembly(args):
    if args.track:
        createtrackfile('assemble',args.nodes)
    it = 1
    while it <= args.nodes:
        count = (it-1)*args.cores+1+args.startcount
        update_hyakrun(args,count)
        os.system('qsub -q bf hyakrun.sh')
        it+=1

#run the assembly
def runassembly(count):
    os.system('sh 1C.sh %s' % count)

#run the assembly in parallel
def runparallelassembly(args):
    jobcount = 1
    jobid = args.startcount
    while jobcount <= args.cores:
          if __name__ == '__main__':
             jobs = []
             for i in range(1, args.cores+1):
                 p = multiprocessing.Process(target=runassembly, args=(jobid,))
                 jobs.append(p)
                 p.start() 
                 jobcount = jobcount + 1
                 jobid +=1
                 if jobcount > args.cores:
                     break
             for x in jobs:
                 x.join()
    if args.track:
        updatetrackfile()

#run a single docking job
def rundock(residue):
    os.system('sh 1A.sh %s' % residue)

#run docking jobs on multiple cores
def runparalleldocking(args):
    residues = args.residues
    compnumber = args.cores
    jobcount = 0
    while jobcount < len(residues):
        if __name__ == '__main__':
            jobs = []
            it = 0
            while it < compnumber:
               if not os.path.isfile('round%s/input.%s.silent'%(args.roundcount,residues[jobcount])):
                    p = multiprocessing.Process(target=rundock, args=(residues[jobcount],))
                    jobs.append(p)
                    p.start()
                    it+=1
               jobcount = jobcount + 1
               if jobcount >= len(residues):
                   break
            for x in jobs:
               x.join()
    if args.track == True:
        updatetrackfile()

#Read the fasta in 1A.sh and check for unfinished jobs. Que as many nodes as necessary to complete.
def setupdocking(args):
    fasta = ''
    os.system('mkdir round%s'%args.roundcount)
    get_fasta_from_script(args)
    totalresidues = count_residues(args.fasta)
    fragsize = checkfragsize()
    totaljobs = totalresidues-fragsize+1
    newjobs = find_missing(args,totaljobs)
    totaljobs = len(newjobs)
    totalnodes = int(totaljobs / args.cores) + (totaljobs % args.cores > 0)
    if args.track == True:
        createtrackfile('dock',totalnodes)
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
        update_hyakrun_dock(args,jobresidues)
        os.system('qsub -q bf hyakrun.sh')
        lower = upper

def get_fasta_from_script(args):
    fasta = ''
    with open('1A.sh', 'r') as runscript:
        for line in runscript:
            if 'in::file::fasta' in line:
                line = line.strip().split()
                fasta = line[1]
    args.fasta

#get the number of residues in the fasta
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

#Check the size of fragments being used
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
                        if 'position' in line:
                            robettastyle = True
                            continue
                        if robettastyle == True and line.strip() == '' and firstspace == True:
                            firstspace = False
                            continue
                        if robettastyle == True and line.strip() == '' and firstspace == False:
                            return fragsize
                        if robettastyle == True and firstspace == False:
                            fragsize +=1
                    return fragsize

#check the round folder for existing silentfiles
def find_missing(args,expectedtotal):
    missing = []
    for i in range(1, expectedtotal+1):
        if os.path.isfile('round%s/input.%d.silent' % (args.roundcount,i)) == False:
            missing.append(i)
    print missing
    return missing

#count the calphas in the pdb. Only counts fragments greater than 9 residues
def count_calphas(pdb):
    cacount=0
    with open(pdb,'r') as pdb:
        currentfragcount = 0
        previousres = 0
        for line in pdb:
            if line.startswith('ATOM') and 'CA' in line:
                resnum = int("".join(list(line)[22:26]))
                if resnum == previousres+1:
                    currentfragcount+=1
                else:
                    if currentfragcount >= 9:
                        cacount+=currentfragcount
                    currentfragcount = 0
                previousres = resnum
        if currentfragcount >= 9:
            cacount+=currentfragcount
    return cacount

#return the percentage of the residues that have been assigned (only counts 9+ residue fragments)
def percent_finished(args):
    fasta = ''
    with open('1A.sh', 'r') as runscript:
       for line in runscript:
           if 'in::file::fasta' in line:
                line = line.strip().split()
                fasta = line[1]
    totalresidues = count_residues(fasta)
    assignedresidues = count_calphas('round%s_model.pdb'%args.roundcount)
    percent = 100*float(assignedresidues)/float(totalresidues)
    return percent

#return the total residues that have been added to the pose since the last round (9+ fragments only)
def total_added(args):
    startcount = 0
    if args.roundcount > 1:
        startcount = count_calphas('round%sinput.pdb'%args.roundcount)
    endcount = count_calphas('round%s_model.pdb'%args.roundcount)
    totaladded = endcount-startcount
    return totaladded

#functions for post docking preparation
def grower_prep(args):
    pdb = open('round%s_model.pdb'%args.roundcount,'r').readlines()
    frags = convert_pdb_to_frags(pdb)
    setup_grower_pdb_and_coords(args,frags)
    total_loops = read_growerprep()
    if not os.path.isdir('grower/runfiles/'):
        print "no grower directory. Exiting protocol"
        exit()
    os.system('cp *.mers grower/runfiles/')
    os.system('cp input.pdb grower/runfiles/')
    os.chdir('grower/runfiles/')
    gcommand = ('python multinodegrower.py -name %s -lc %s'%(args.name,total_loops))
    os.system(gcommand)


def setup_grower_pdb_and_coords(args,frags):
    for frag in frags:
        fsize = fragCAcount(frag)
        if fsize < args.fragsize:
            add_frag_to_coordfile(frag)
        else:
            write_frag_to_pdb(args,frag)

def read_growerprep():
    with open('grower_prep.txt','r') as gp:
        for line in gp:
            if 'total loops' in line:
                line = line.split()
                total_loops = int(line[3])
                return total_loops

#read the pdb and split it into fragments
def convert_pdb_to_frags(pdb):
    pdbfrags = []
    atomlist = []
    previousresnum = 0
    resN = []
    resNsub1 = []
    previousC = ''
    currentN = ''
    firstres = True
    firstatom = True
    for line in pdb:
        if line.startswith('ATOM'):
            line = list(line)
            atomid = "".join(line[11:16]).strip()
            resid = line[17:20]
            resnum = int("".join(line[22:26]))
            x = float("".join(line[31:39]))
            y = float("".join(line[39:47]))
            z = float("".join(line[47:54]))
            atom = Atom("".join(line),resnum,atomid,resid,x,y,z)
            #handle the first atom/residue specifically
            if firstatom == True:
                resN.append(atom)
                firstatom = False
            elif firstres == True:
                if previousresnum != resnum:
                    resNsub1 = [x for x in resN]
                    for atomit in resNsub1:
                        atomlist.append(atomit)
                    del resN[:]
                    firstres = False
                resN.append(atom)
            else:
                #on new res check distance and shift one res in the storage.
                if previousresnum != resnum:
                    for atomiter in resN:
                        if atomiter.atomid == 'N':
                            currentN = atomiter
                    for atomiter2 in resNsub1:
                        if atomiter2.atomid == 'C':
                            previousC = atomiter2
                    resNsub1 = [x for x in resN]
                    resN = []
                    for atomiter in resNsub1:
                        atomlist.append(atomiter)
                if resnum-previousresnum > 1:
                    fragment = Fragment(atomlist)
                    del atomlist[:]
                    pdbfrags.append(fragment)
                resN.append(atom)
            previousresnum = resnum
    for atomiter in resN:
        atomlist.append(atomiter)
    fragment = Fragment(atomlist)
    pdbfrags.append(fragment)
    
    return pdbfrags

def add_frag_to_coordfile(frag):
    with open('coordfile.txt','a') as coordfile:
        for atom in frag.atoms:
            if atom.atomid == "CA":
                coordline = str(atom.resnum) + " 2 1 " + str(atom.x) + " " + str(atom.y) + " " + str(atom.z)+"\n"
                coordfile.write(coordline)

def write_frag_to_pdb(args,frag):
    with open('input.pdb','a') as pdbfile:
        totalatoms = fragCAcount(frag)
        count = 0
        currentatoms = []
        for atom in frag.atoms:
            if atom.atomid == 'N':
                count+=1
            if count > args.trim and count <= totalatoms-args.trim:
                pdbfile.write(atom.line)

def fragCAcount(frag):
    rescount = 0
    for atom in frag.atoms:
        if atom.atomid == 'CA':
            rescount+=1
    return rescount

def que_fragpicking(args):
   update_hyakrun(args,0)
   os.system('qsub -q bf hyakrun.sh')

def run_fragpicking(args):
    os.system('sh growerprep.sh input.pdb %s'%args.fasta)
    if args.track:
        updatetrackfile()

class Atom:
    def __init__(self,line,num,atomid,code,x,y,z):
        self.line = line
        self.resnum = num
        self.atomid = atomid
        self.code = code
        self.x = x
        self.y = y
        self.z = z

class Fragment:
    def __init__(self,atoms):
        self.atoms = [x for x in atoms]

main()
