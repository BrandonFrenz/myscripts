#!/usr/bin/python
import argparse
from amino_acids import longer_names 

def main():
    args = parseargs()
    if args.rep_end - args.rep_start != args.pos_end-args.pos_start:
        print " ranges don't match"
        exit()
    coordpdb = open(args.positionpdb,'r').readlines()
    noncoordpdb = open(args.pdb,'r').readlines()
    atoms = grabatoms(args,coordpdb)
    movedpdb = copycoords(args,noncoordpdb,atoms)
    print "".join(movedpdb)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',help='The pdb you want to use as a base')
    parser.add_argument('-pp','--positionpdb',help='The pdb you want to steal the coordinates from')
    parser.add_argument('-pps','--pos_start',type=int,help='The start number of the coords you want to steal')
    parser.add_argument('-ppe','--pos_end',type=int,help='The end number of the coords you want to steal')
    parser.add_argument('-rs','--rep_start',type=int,help='The start residue of the positions you want to replace')
    parser.add_argument('-re','--rep_end',type=int,help='The end residue of the position you want to replace')
    args = parser.parse_args()
    return args

def get_amino_acid_lines(pdb):
    atomlines = []
    for line in pdb:
        if line.startswith('ATOM'):
            ll = list(line)
            resid = ''.join(ll[17:20])
            if resid.strip() in longer_names: 
                atomlines.append(line)
    return atomlines

def get_sugar_atom_lines(pdb):
    suglines = []
    for line in pdb:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            ll = list(line)
            resid = ''.join(ll[17:20])
            if resid.strip() not in longer_names: 
                suglines.append(line)
    return suglines

def get_hetnam_lines(pdb):
    nonatomlines = []
    for line in pdb:
        if line.startswith('HETNAM'):
            nonatomlines.append(line)
    return nonatomlines

def get_link_lines(pdb):
    linklines = []
    for line in pdb:
        if line.startswith('LINK'):
            linklines.append(line)
    return linklines

def get_trailing_lines(pdb):
    trail_lines = []
    trailing = False
    for line in pdb:
        if line.startswith('CONECT'):
            trailing = True
        if trailing:
            trail_lines.append(line)
    return trail_lines

#grabs coordinates all all atoms and labels them using rosetta numbering 
def grabatoms(args,coordpdb):
    atomlist = []
    for line in coordpdb:
        if line.startswith('ATOM'):
            ll = list(line)
            atomid = "".join(ll[11:16]).strip()
            resnum = int("".join(line[22:26]))
            if resnum >= args.pos_start and resnum <= args.pos_end:
                resid = ll[17:20]
                x = ll[31:39]
                y = ll[39:47]
                z = ll[47:54]
                atom = Atom(line,resnum,atomid,resid,x,y,z)
                atomlist.append(atom)
    return atomlist

def copycoords(args,pdb,atoms):
    newpdb = []
    atomcount = 0
    for line in pdb:
        if line.startswith('ATOM'):
            ll = list(line)
            atomid = "".join(ll[11:16]).strip()
            current_resnum = int("".join(ll[22:26]))
            if current_resnum >= args.rep_start and current_resnum <= args.rep_end:
                atom = atoms[atomcount]
                atomcount +=1
                ll[31:39] = atom.x
                ll[39:47] = atom.y
                ll[47:54] = atom.z
            line = "".join(ll)
        newpdb.append(line)
    return newpdb

class Atom:
    def __init__(self,line,num,atomid,code,x,y,z):
        self.line = line
        self.resnum = num
        self.atomid = atomid
        self.code = code
        self.x = x
        self.y = y
        self.z = z

main()
