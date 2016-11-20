#!/usr/bin/python
import argparse
import re
import pdbtools

def main():
    args = parseargs()
    args.residues = parse_residues(args)
    remove_residues(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r','--residues',nargs='+',help='The residues to remove')
    parser.add_argument('-p','--pdbs',nargs='+',help='The pdb(s) to remove the residues from')
    parser.add_argument('-o','--output',help='The name of the output, if more than 1 pdb outputs will have an iterator applied. To overwrite use "0"')
    args = parser.parse_args()
    return args

def parse_residues(args):
    residues = []
    for res in args.residues:
        res = res.split('-')
        chain = re.split('\d+',res[0])[0]
        lower = int(re.split('(\d+)',res[0])[1])
        upper = lower
        if len(res) > 1:
            upper = int(res[1])
        for i in range(lower,upper+1):
                residues.append(i)
    return residues

def remove_residues(args):
    it = 1
    for pdb in args.pdbs:
        pdbfile = open(pdb,'r').readlines()
        residues = pdbtools.get_residue_list(pdbfile)
        newres = []
        for res in residues:
            if res.num in args.residues:
                continue
            else:
                newres.append(res)

        residues = newres
        newpdb = pdbtools.make_pdblines_from_residues(residues)
        nameit = it
        if len(args.pdbs) == 1:
            nameit = ''
        pdbname = args.output+"%s"%nameit
        if args.output == '0':
            pdbname = pdb
        pdbtools.write_pdb(newpdb,pdbname)

main()
