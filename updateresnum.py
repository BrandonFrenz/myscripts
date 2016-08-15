#!/usr/bin/python
import argparse
import re

def main():
    args = parseargs()
    for pdb in args.pdbs:
        newpdb = updatepdb(args,pdb)
        writepdb(args, pdb,newpdb)


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdbs', nargs='+', help='the pdbs structures to update')
    parser.add_argument('-a','--add', type=int, help='the number to increase each residue by')
    parser.add_argument('-o','--overwrite', action='store_true', help='overwrite the old pdb name with the updated pdbfile')
    args = parser.parse_args()
    return args

def updatepdb(args,pdbfile):
    newpdb = []
    with open(pdbfile,'r') as pdb:
        for line in pdb:
            if line.startswith('ATOM'):
                residue = line[22:26]
                rescount = int(residue.strip())
                newrescount = rescount+args.add
                newressub = " "*(4-len(str(newrescount)))+str(newrescount)
                listline = list(line)
                listline[22:26]=newressub
                line = "".join(listline)
            newpdb.append(line)
    return newpdb  

def writepdb(args, pdb, newpdb):
    name = "updated_"+pdb
    if args.overwrite:
        name = pdb
    with open(name,'w') as newpdbfile:
        for line in newpdb:
            newpdbfile.write(line)

main()
