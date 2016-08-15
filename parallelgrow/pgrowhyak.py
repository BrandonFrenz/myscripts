#!/usr/local/bin/python2.7
import sys
import argparse
from multiprocessing import Pool
from os import popen,  system
from os.path import basename, exists
from sys import exit, stderr, stdout
import operator
import os
import multiprocessing
import time
import glob
import fileinput
import re
import shutil

def main():
    args = parseargs()
    if args.mode == 'prep':
        setup_and_que(args)
    if args.mode == 'run':
        if not os.path.isfile('finished.txt'):
            growloop(args) #grow all the loops each loop needs a directory with it's number in the sequence alignment
        postloopstart = time.time()
        if 'context' not in os.getcwd() and args.compare:
            args.context = False
        else:
            args.context = True
        make_lps_from_taboo(args)
        os.chdir('../')
        multiplerounds = True ####this generates new lps files in that includes the results of all the previous runs if there is only 1 run it's a waste of time and you can turn it off
        haslps = check_lps(args) # Moves up a directory and checks to see that each loop has finished
        print haslps
        if haslps:
            collect_lps(args) #combines the loop partial solution files for each loop and moves them to the comparator directory if you've created it ahead of time
            runcomparator(args) #If a directory named 'comparator' exists it will attempt to combine the loops. You need to create the 'comparator' directory manually if you want this step'
            needsanother = needsmorerounds() #Checks if the models clash and if not the next steps will relax the model and if a native structure is provided compare the RMSDs to the native.
            if needsanother == False:
                relaxpdbs(args)
                checkrms(args)
                prepcontext(args)
        postlooptime = time.time()-postloopstart
        os.system('echo "%s" >> time.txt'%postlooptime)
    if args.mode == 'prepcontext':
        prepcontext(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m','--mode',default='prep',help='The mode prep or run')
    parser.add_argument('-p','--pdb',default='input.pdb', help='The input pdb')
    parser.add_argument('-x','--xml',default='parallel.xml',help='The xml to run')
    parser.add_argument('-d','--density',default='density.mrc',help='The electron density map')
    parser.add_argument('-n','--native',default='native.pdb',help='The native pdb')
    parser.add_argument('-name','--name',help='The name of the job')
    parser.add_argument('-l','--loop',type=int,default=0,help='Which loop are you operating on. Default will leave xml unchanged')
    parser.add_argument('-c','--cores',type=int,default=16,help='How many cores you have available per node')
    parser.add_argument('-wt','--walltime',type=int,default=5,help='What do you want the wall time to be?')
    parser.add_argument('-t','--taboo',type=int,default=1,help='Do you want to use taboo filtering when available')
    parser.add_argument('-compare','--compare',type=bool,default=True,help='Do you want to run the comparator when all jobs are finished?')
    parser.add_argument('-context','--context',type=bool,default=False,help='Is this context growing or not? This will be set based on directory name')
    args = parser.parse_args()
    return args

#enter the directory args.name setup a job and que it
def setup_and_que(args):
    if not os.path.isdir(args.name):
        os.mkdir(args.name)
    for runfile in glob.glob('*'):
        if not os.path.isdir(runfile):
            shutil.copy(runfile,args.name)
    os.chdir(args.name)
    update_hyakscript(args)
    os.system('qsub -q bf hyakrun.sh')

def update_hyakscript(args):
    cwd = os.getcwd()
    for line in fileinput.input('hyakrun.sh', inplace=1):
       if line.startswith('#PBS'):
            line = re.sub("-N \S+", "-N %s" % (args.name), line)
            line = re.sub("walltime=\S+", "walltime=%s:00:00"%(args.walltime), line)
            line = re.sub("ppn=\d+,", "ppn=%s,"%(args.cores), line)
            line = re.sub("feature=\d+", "feature=%s"%(args.cores), line)
            line = re.sub("-d \S+", "-d %s/" % (cwd), line)
       if 'pgrowhyak.py' in line:
            line = "python pgrowhyak.py -m run -x %s -p %s -d %s -n %s -c %s -t %s &> /dev/null \n" % (args.xml, args.pdb, args.density, args.native, args.cores, args.taboo)
       sys.stdout.write(line)
    
    if args.loop != 0:
        for line in fileinput.input('parallel.xml', inplace=1):
           sys.stdout.write(re.sub("looporder=\S+", "looporder=%s"%args.loop, line))

def growloop(args):
    if os.path.isfile('finished.txt'):
        print 'found finished'
        return
    if args.taboo == 1:
        if not os.path.isdir('../taboo'):
            os.mkdir('../taboo')
    counter = get_counter()

    
    start_time = time.time()
    #this first run generates doesn't need to be parallel and it generates the files required for the first round
    cleanprevious()
    
    #this parameter determins whether initial setup needs to be done on the first pass. It doesn't if the job isn't starting from the beginning
    ready_for_parallel = False
    if counter == 1:
        #I'm currentlly skipping taboo on first beam.
        filterbeam = 0;
        beamfile = 'none'
        #if os.path.isfile('beam1.txt'):
        #    filterbeam = 1
        #    beamfile = 'beam1.txt'
        os.system("sh prun.sh %s %s %s %s 0 NA 1 1 %s %s" % (args.xml,args.pdb,args.density,args.native,filterbeam,beamfile) )
    else:
        splitjobs(counter,args.cores)
        ready_for_parallel = True



    done = False
    if os.path.isfile('finished.txt'):
        print 'found finished'
        done = True
    setuptime = time.time()-start_time
    os.system('echo %s >> time.txt'%setuptime)
    while done == False:
        begintime = time.time()
        if ready_for_parallel == False:
            combine(counter)
            beamcount = secondfilter(args,counter)
            splitjobs(counter,args.cores)
        ready_for_parallel = False
        runparalleljobs(args,counter)
        if args.taboo == 1:
            addtaboo(counter)
        counter+=1
        count = open('count', "w")
        count.write('%d' % (counter),)
        count.close()
        if(os.path.isfile("./finished.txt")):
           done = True
           print "found finished"
        roundtime = (time.time()-begintime)*args.cores#this won't be totally accurate if you don't actually use all the cores.
        os.system('echo %s >> time.txt'%roundtime)
    #this job is used to print the lps files at the end
    finalstepstarttime = time.time()
    combine(counter)
    beamcount = secondfilter(args,counter)
    if args.taboo == 1:
        addtaboo(counter)
   
    totaltime = time.time()-start_time
    finalsteptime = time.time()-finalstepstarttime
    os.system('echo %s >> time.txt'%finalsteptime)
    print "the total time to run is "
    print totaltime

#this file gets the current beam count by reading the beam_0.txt file (it no longer reads the count file which isn't necessary)
def get_counter():
    counter = 1
    if os.path.exists('beam_0.txt'):
       with open('beam_0.txt') as f:
           counter = int(f.readline().split()[-2])
    return counter

#remove all beams
def cleanprevious():
    for jobfiles in glob.glob("beam_*_*.txt"):
        os.remove(jobfiles)
    for jobfiles in glob.glob('beam_*.*.txt'):
        os.remove(jobfiles)

def run(args,filename, jobcounter):
    previouscount = int(filename.split('_')[1])+1
    previousbeams = '../taboo/beam%s.txt'%(previouscount)
    filterbeam = 0
    if args.taboo == 1:
        if checktaboo(previouscount):
            filterbeam = 1
    command = "sh prun.sh %s %s %s %s 1 %s 1 %s %s %s" % (args.xml,args.pdb,args.density,args.native,filename,jobcounter, filterbeam, previousbeams) 
    os.system(command)
    print " just ran  " + command

#merge the results from each core
def combine(counter):
    combinedbeams = open("beam_%d.txt" % (counter), "w")
    for beamfile in glob.glob('beam_%d.*.txt'%counter):
        singleroundbeam = open(beamfile, "r")
        combinedbeams.write(singleroundbeam.read())
        singleroundbeam.close()
    combinedbeams.close()
    cleanprevious()

#filter the combined results
def secondfilter(args,counter): 
    filename = ("beam_%d.txt" % (counter))
    filterbeam = 0
    beamfile = 'none'
    if args.taboo == 1:
        if checktaboo(counter):
            filterbeam = 1
            beamfile = '../taboo/beam'+str(counter)+'.txt'
    os.system("sh prun.sh %s %s %s %s 1 %s 1337 0 %s %s" % (args.xml,args.pdb,args.density,args.native,filename, filterbeam, beamfile) )
    print "finished sorting " + filename
    storedbeam = open("beam_0.txt", "r")
    beamcount = 0
    for line in storedbeam.readlines():
        if len(line.split()) == 5: 
            beamcount+=1
    storedbeam.close()
    print beamcount
    return beamcount

#split the beamfiles up to be run on all the cores
def splitjobs(counter,cores):
    storedbeam = open("beam_0.txt", "r")
    divisioncount = 1
    partialset = open("beam_0.txt", "r")
    newbeam_linesize = 5
    for line in storedbeam.readlines():
        if len(line.split()) == newbeam_linesize:
            partialset.close()
            partialset = open("beam_%d_%d.txt" % (counter,divisioncount), "a")
            divisioncount+=1
        if divisioncount > cores:
            divisioncount = 1
        partialset.write(line)
    storedbeam.close()    
    partialset.close()   

#run the loop grower in parallel
def runparalleljobs(args,counter):
    jobcount = 1
    filecount = len(glob.glob('beam_%d_*.txt'%counter))
    while jobcount <= filecount:
          if __name__ == '__main__':
             jobs = []
             for i in range(1, args.cores+1):
                 p = multiprocessing.Process(target=run, args=(args,"beam_%d_%d.txt" % (counter,i),i))
                 jobs.append(p)
                 p.start() 
                 jobcount = jobcount + 1
                 if jobcount > filecount:
                     break
             for x in jobs:
                 x.join()

def checktaboo(count):
    return os.path.isfile('../taboo/beam%s.txt'%count)

#Check all the beams in the currentfile and in the taboo file. Add all unique beams to the taboo file.
def addtaboo(count):
    beamlist = []
    firstline = ''
    with open('beam_0.txt','r') as beamfile:
        beam = []
        for line in beamfile:
            if firstline == '':
                firstline = line
            if line == firstline:
                if len(beam) != 0:
                    beamlist.append(beam)
                    beam = []
            beam.append(line)
        beamlist.append(beam)
    
    taboobeamlist = []
    firstline = ''
    if checktaboo(count):
        with open('../taboo/beam%s.txt'%count) as taboobeams:
            beam = []
            for line in taboobeams:
                if firstline == '':
                    firstline = line
                if line == firstline:
                    if len(beam) != 0:
                        taboobeamlist.append(beam)
                        beam = []
                beam.append(line)
            taboobeamlist.append(beam)
    uniquebeams = []
    for beam in beamlist:
        if beam not in taboobeamlist:
            uniquebeams.append(beam)
    with open('for_taboo.txt','w') as fortaboo:
        for beam in uniquebeams:
            for line in beam:
               fortaboo.write(line)
    os.system('cat for_taboo.txt >> ../taboo/beam%s.txt'%count)

def count_beams(beamfile):
    firstline = ''
    beamcount = 0
    with open(beamfile,'r') as bfile:
        for line in bfile:
            if firstline == '':
                firstline = line
            if firstline == line:
                beamcount+=1
    return beamcount

def check_lps(args):
    havelpsfiles = True
    loopdirs = []
    if args.context == True:
        loopdirs = glob.glob('context[0-9]*')
    else:
        loopdirs = glob.glob('[0-9]*')
    for i in loopdirs:
       loopnum = i
       if args.context == True:
           loopnum = 1
       print i
       lpsfile = i+"/lpsfile_%s.0.txt"%loopnum
       if not os.path.isfile(lpsfile):
           havelpsfiles = False

    return havelpsfiles


#this function generates lps files from the stored beams. If you don't have stored beams running it won't hurt other than wasting a few minutes of cpu time.
def new_lps_file(args):
    #update xml
    for line in fileinput.input(args.xml, inplace=1):
           sys.stdout.write(re.sub("beamwidth=\d+", "beamwidth=250", line))
    
    maxcount = 0
    for i in glob.glob('beam[0-9]*.txt'):
        split_i = re.split('beam|.txt',i)
        count = int(split_i[1])
        if count > maxcount:
            maxcount = count
    command = "sh prun.sh %s %s %s %s 1 beam%s.txt 1337 0 0 na"%(args.xml,args.pdb,args.density,args.native,maxcount)
    os.system(command)


def collect_lps(args):
    with open("lpsfile.txt", "w") as combinedlps:
        loopcount = 0
        loopdirs = glob.glob('[0-9]*')
        if args.context == True:
            loopdirs = glob.glob('context[0-9]*')
        stringloopcount = str(len(loopdirs))
        combinedlps.write(stringloopcount+ "\n")
        for loopdir in loopdirs:
            beamcount = 0;
            for x in glob.glob('%s/lpsfile*' % (loopdir)):
                lpsfile = open(x,"r")
                for line in lpsfile.readlines():
                    if len(line.split()) == 1:
                        beamcount = beamcount + 1
                beamcount = beamcount/2        
                strbeamcount = str(beamcount)
                combinedlps.write(strbeamcount + "\n")
                lpsfile.seek(0)
                combinedlps.write(lpsfile.read())
    if args.context == True:
        if os.path.isdir('contextcomparator'):
            shutil.copy('lpsfile.txt','contextcomparator')
    else:
        if os.path.isdir('comparator'):
            shutil.copy('lpsfile.txt','comparator')

def runcomparator(args):
    compdir = 'comparator'
    if args.context == True:
        compdir = 'contextcomparator'
    if os.path.isdir(compdir):
       os.chdir(compdir)
       if os.path.isfile('startedcomparator'):
           exit()
       runningcomparator = open('startedcomparator','w')
       runningcomparator.close()
       for file in glob.glob(os.path.join('../runfiles/', '*')):
         shutil.copy(file, './')
       for line in fileinput.input(args.xml, inplace=1):
           sys.stdout.write(re.sub("read_from_file=0", "read_from_file=1", line))
       command = "sh prun.sh %s %s %s %s 0 na 0 0 0 na" % (args.xml,args.pdb,args.density,args.native)
       os.system(command)

def checkrms(args):
    if os.path.isfile('native.pdb'):
        jobcount = 0
        pdblist = glob.glob('aftercomparator*0001.pdb')
        pdbcount = len(pdblist)
        while jobcount < pdbcount:
            if __name__ == '__main__':
                jobs = []
                for i in range(0,args.cores):
                    pdb = pdblist[jobcount]
                    p = multiprocessing.Process(target=runrms, args=("%s" % (pdb),) )
                    jobs.append(p)
                    p.start() 
                    jobcount = jobcount + 1
                    if jobcount >= pdbcount:
                        break
                for i in jobs:
                    i.join()

def runrms(pdb):
    command = ('/gscratch/dimaio/bfrenz/Rosetta/score_jd2.static.linuxgccrelease -in:file:native native.pdb -s %s -ignore_unrecognized_res -edensity:mapfile density.mrc'
    ' -fastdens_wt 20 -out:suffix rms'%pdb)
    os.system('%s' % command)

def needsmorerounds():
    needs_another = True
    if os.path.isfile('recommendation.txt'):
        with open('recommendation.txt') as recommendation:
            for line in recommendation:
                if '''Another round probably isn't necessary''' in line:
                    needs_another = False
    return needs_another

def runrelax(pdb):
    command = ('/gscratch/dimaio/bfrenz/Rosetta/rosetta_scripts.static.linuxgccrelease -s %s -parser:protocol relax.xml -edensity:mapfile density.mrc -default_max_cycles 200'
        ' -missing_density_to_jump -beta_july15'%pdb)
    os.system('%s' % command)

def relaxpdbs(args):
    if os.path.isfile('relax.xml'):
        jobcount = 0
        pdblist = glob.glob('aftercomparator*')
        pdbcount = len(pdblist)
        while jobcount < pdbcount:
            if __name__ == '__main__':
                jobs = []
                for i in range(0,args.cores):
                    pdb = pdblist[jobcount]
                    p = multiprocessing.Process(target=runrelax, args=("%s" % (pdb),) )
                    jobs.append(p)
                    p.start() 
                    jobcount = jobcount + 1
                    if jobcount >= pdbcount:
                        break
                for i in jobs:
                    i.join()
    else:
        print "you are missing the xml for relax"
        exit()

def prepcontext(args):
    createcontextposes()
    maxcontext = 0
    for i in glob.glob('context*'):
        count = re.split('context|.pdb',i)[1]
        if int(count) > maxcontext:
            maxcontext = int(count)
    print maxcontext

    count = 1
    os.chdir('../')
    while count <= maxcontext:
        contextdir = 'context'+str(count)
        if os.path.isdir(contextdir):
            count+=1
            continue
        else:
            os.mkdir(contextdir)
            os.chdir(contextdir)
            os.system('cp ../runfiles/* ./')
            os.system('cp ../comparator/context%s.pdb ./'%(count))
            args.loop = 1
            update_hyakscript(args)
            os.chdir('../')
        count+=1

def readloopfile():
    loops = {}
    with open ('loopresidues.txt','r') as loopresidues:
        for line in loopresidues:
            line = line.split()
            loopnum = line[1]
            loopstart = line[2]
            loopstop = line[3]
            loops[loopnum] = [loopstart,loopstop]
    return loops

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def lowestscore():
    bestscore = 9999999
    bestname = ''
    with open('score.sc','r') as scores:
        for line in scores:
            if line.startswith('SCORE:'):
                line = line.split()
                if isfloat(line[1]):
                    score = float(line[1])
                    name = line[len(line)-1]
                    if float(score) < bestscore:
                        bestscore = score
                        bestname = name
    
    bestname = bestname+".pdb"
    return bestname

def cutpose(name,loopnum,start,stop):
    newpdb = []
    with open(name,'r') as pdb:
        for line in pdb:
            if line.startswith('ATOM'):
                res = line[23:26]
                res = int(res)
                if res < start or res > stop:
                    newpdb.append(line)
    with open('context%s.pdb'%loopnum,'w') as pdbfile:
        for i in newpdb:
            pdbfile.write(i)
            

def createcontextposes():
    loops = readloopfile()
    lowestname = lowestscore()
    for loop in loops.keys():
        start = int(loops[loop][0])
        stop = int(loops[loop][1])
        cutpose(lowestname,loop,start,stop)

def make_lps_from_taboo(args):
    for line in fileinput.input(args.xml, inplace=1):
           sys.stdout.write(re.sub("beamwidth=\d+", "beamwidth=250", line))
    os.system('cp parallel.xml ../')
    os.chdir('../')
    beams = glob.glob('taboo/beam*')
    maxbeam = 0
    for beam in beams:
        beam = int(re.split('beam|.txt',beam)[1])
        if beam > maxbeam:
            maxbeam = beam
    #print maxbeam
    command = ('sh prun.sh %s %s %s %s 1 taboo/beam%d.txt 1337 0 0 na')%(args.xml,args.pdb,args.density,args.native,maxbeam)
    os.system(command)
    os.system('rm beam_0.txt')

main()
