#!/usr/bin/python
import argparse
import re
import numpy
import math
import os

def main():
    args = parseargs()
    inputmap = mapfrags(args.pdbs)
    if args.mode == 'prep':
        inputmap = prepareinputs(inputmap,args.sizecut)
        writeall(inputmap)
    if args.mode == 'stats':
        inputmap = prepareinputs(inputmap,args.sizecut)
        inputmap = allfragrms(args,inputmap)
        reportstats(inputmap)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m','--mode',help='Which functions to call')
    parser.add_argument('-p','--pdbs',nargs="+",help='The pdbs to processes')
    parser.add_argument('-n','--natives',nargs="+",help='The native structures')
    parser.add_argument('-f','--fastas',nargs="+",help='The fastas')
    parser.add_argument('-sc','--sizecut',type=int,default=15,help='Cut all fragments shorter than this')
    args = parser.parse_args()
    return args

#split the pdb into a list of fragment objects
def readpdb(pdb):
    pdbfrags = []
    with open(pdb,'r') as pdbfile:
        atomlist = []
        previousresnum = 0
        resN = []
        resNsub1 = []
        previousC = ''
        currentN = ''
        firstres = True
        firstatom = True
        for line in pdbfile:
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
                    dist = 0
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
                        dist = distancecheck(previousC,currentN)
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

def count_CApdb(pdb):
    cacount = 0
    with open(pdb,'r') as pdbfile:
        for line in pdbfile:
            if line.startswith('ATOM'):
                if 'CA' in line:
                    cacount+=1
    return cacount

def count_fasta(fasta):
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

def verifynatives(args):
    for native in args.natives:
        pdbcount = count_CApdb(native)
        key = re.split('/|native.pdb',native)[-2]
        for fasta in args.fastas:
            if key in fasta:
                fastacount = count_fasta(fasta)
                if fastacount != pdbcount:
                    print "Warning " + native + " does not have the same number of residues as the fasta"
                    print "native has " + str(pdbcount) + " fasta has " + str(fastacount)

def mapfrags(pdblist):
    fragmap = {}
    for pdb in pdblist:
        frags = readpdb(pdb)
        fragmap[pdb] = frags
    return fragmap

def checkbreaks(args):
    for native in args.natives:
        frags = readpdb(native)
        if len(frags) > 1:
            print "breaks in " + str(native)
            chainbreaks = []
            isfirst = True
            for frag in frags:
                if not isfirst:
                    chainbreaks.append(str(frag.atoms[-1].resnum))
                isfirst = False
            print " ".join(chainbreaks)

                
def distancecheck(atom1,atom2):
    a1 = numpy.array((atom1.x,atom1.y,atom1.z))
    a2 = numpy.array((atom2.x,atom2.y,atom2.z))
    dist = numpy.linalg.norm(a1-a2)
    return dist

def fragrms(frag,native):
    totaldist = 0
    for atom in frag.atoms:
        if atom.atomid != 'CA':
            continue
        natatom = get_CA(native,atom.resnum)
        if natatom is None:
            return -1
        dist = distancecheck(atom,natatom)
        totaldist += dist*dist
    rms = math.sqrt(totaldist/len(frag.atoms))
    return rms

def get_CA(pdb,targetres):
    with open(pdb,'r') as pdbfile:
        realrescount = 0
        for line in pdbfile:
            if line.startswith('ATOM'):
                line = list(line)
                atomid = "".join(line[11:16]).strip()
                if atomid != "CA":
                    continue
                realrescount+=1
                resid = line[17:20]
                resnum = int("".join(line[22:26]))
                if realrescount != targetres:
                    continue
                x = float("".join(line[31:39]))
                y = float("".join(line[39:47]))
                z = float("".join(line[47:54]))
                atom = Atom("".join(line),resnum,atomid,resid,x,y,z)
                return atom

def allfragrms(args,pdbmap):
    for pdb in pdbmap.keys():
        pdbkey = re.split('_',pdb)[0]
        for native in args.natives:
            #nativekey = re.sub('-','',native)
            if pdbkey in native:
                it = 0
                while it < len(pdbmap[pdb]):
                    frag = pdbmap[pdb][it]
                    frms = fragrms(frag,native)
                    pdbmap[pdb][it].rms = frms
                    it+=1
    return pdbmap

def reportstats(pdbmap):
    longestbadfrag = ''
    longestbadfraglen = 0
    print ' '
    for pdb in sorted(pdbmap.keys()):
        for fragment in pdbmap[pdb]:
            rms = -1
            if hasattr(fragment,'rms'):
                rms = fragment.rms
            fragreport =  str(pdb) + " " + str(fragment.atoms[0].resnum) + "-" + str(fragment.atoms[-1].resnum) + " " + str(rms)
            print fragreport
            if rms > 2.0:
                if len(fragment.atoms) > longestbadfraglen:
                    longestbadfraglen = len(fragment.atoms)
                    longestbadfrag = fragreport
        print ' '
    if longestbadfrag == '':
        print " Congratulations there are no bad fragments "
    else:
        print "the longest bad fragment is " + str(longestbadfraglen) + " the offender is:"
        print longestbadfrag
                    
def cleanandmovenatives(args):
    for native in args.natives:
        nativekey = re.sub('-','',native)
        nativekey = nativekey.split('/')[1]+'native.pdb'
        newnative = []
        with open(native,'r') as nativefile:
            for line in nativefile:
                if line.startswith('ATOM'):
                    newnative.append(line)
        with open(nativekey,'w') as newnativefile:
            for line in newnative:
                newnativefile.write(line)

def cleanandmovefastas(args):
    for fasta in args.fastas:
        fastakey = re.sub('-','',fasta)
        fastakey = fastakey.split('/')[1]+'.fasta'
        command = 'cp %s ./%s'%(fasta,fastakey)
        os.system(command)

def writepdbfromfrags(fraglist,name):
    newpdb = []
    lfraglen = 0
    for frag in fraglist:
        for atom in frag.atoms:
            newpdb.append(atom.line)
    with open(name,'w') as newpdbfile:
        for line in newpdb:
            newpdbfile.write(line)


def prepareinputs(inputmap,sizecut):
    for pdb in inputmap.keys():
        frags = inputmap[pdb]
        keptfrags = []
        maxfraglen = 0
        for frag in frags:
            rescount = fragCAcount(frag)
            if rescount >= sizecut:
                keptfrags.append(frag)
            if rescount > maxfraglen:
                maxfraglen = rescount
        if len(keptfrags) == 0:
            for frag in frags:
                rescount = fragCAcount(frag)
                if rescount == maxfraglen:
                    keptfrags.append(frag)
        inputmap[pdb] = keptfrags
    return inputmap
            
def writeall(inputmap):
    for pdb in inputmap.keys():
        frags = inputmap[pdb]
        pdbkey = re.split('_',pdb)[0]
        name = pdbkey+'_nrinput.pdb'
        writepdbfromfrags(frags,name)

def fragCAcount(frag):
    rescount = 0
    for atom in frag.atoms:
        if atom.atomid == 'CA':
            rescount+=1
    return rescount

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

