#!/usr/bin/python
import argparse

def main():
    args = parseargs()
    calphas = grabCAs(args)
    writecalphas(args,calphas)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdbs',nargs='+',help='list of pdbs to grab CAs from')
    parser.add_argument('-n','--name',default='combinedCAs.pdb',help='the name of the CA file')
    args = parser.parse_args()
    return args

def grabCAs(args):
    calphas = []
    for pdb in args.pdbs:
        with open(pdb,'r') as pdbfile:
            for line in pdbfile:
                if line.startswith('ATOM'):
                    atom = line[12:15].strip()
                    if atom == 'CA':
                        calphas.append(line)
    return calphas

def writecalphas(args,calphas):
    with open(args.name,'w') as CAfile:
        for line in calphas:
            CAfile.write(line)

main()
