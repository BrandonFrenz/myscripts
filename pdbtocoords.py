#!/usr/bin/python
import argparse
from collections import defaultdict

def main():
    args = parseargs()
    calphas = get_cas(args)
    write_coords(args,calphas)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdbs',nargs="+",help='The pdbs to utiilize')
    parser.add_argument('-sw','--scoreweight',default='1',help='The score weight to use in the loop growing protocol')
    args = parser.parse_args()
    return args

def get_cas(args):
    calphas = defaultdict(list)
    for pdb in args.pdbs:
        with open(pdb,'r') as pdbfile:
            for line in pdbfile:
                if line.startswith('ATOM'):
                    line = list(line)
                    atomid = "".join(line[11:16]).strip()
                    if atomid != "CA":
                        continue
                    resid = line[17:20]
                    resnum = int("".join(line[22:26]))
                    x = float("".join(line[31:39]))
                    y = float("".join(line[39:47]))
                    z = float("".join(line[47:54]))
                    atom = Atom("".join(line),resnum,atomid,resid,x,y,z)
                    calphas[resnum].append(atom)
    return calphas

def write_coords(args,calphas):
    with open('coordfile.txt','w') as coordfile:
        coordfile.write(args.scoreweight+'\n')
        for res in calphas.keys():
            totalcas = len(calphas[res])
            atomid = "%s 2 "%(res)
            newline = atomid + str(totalcas)
            for ca in calphas[res]:
                newline += " " + str(ca.x) + " " + str(ca.y) + " " + str(ca.z)
            newline += "\n"
            coordfile.write(newline)

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
