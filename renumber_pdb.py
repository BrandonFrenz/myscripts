#!/usr/bin/python
import argparse
import pdbtools

def main():
    args = parseargs()
    renumber_pdb(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',help='The pdbfile')
    parser.add_argument('-o','--output',help='The output file')
    args = parser.parse_args()
    return args

def renumber_pdb(args):
    pdbfile = open(args.pdb).readlines()
    residues = pdbtools.get_residue_list(pdbfile)
    
    #get the header
    headlines = []
    header = True
    for line in pdbfile:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            header = False
        if header == True:
            headlines.append(line)

    tailLines = []
    tailer = True
    it = len(pdbfile)-1
    while it > 0:
        line = pdbfile[it]
        if line.startswith('ATOM') or line.startswith('HETATM'):
            tailer = False
        if tailer == True:
            tailLines.append(line)
        it-=1

    count = 1
    for res in residues:
        res.num = count
        count+=1
    reslines = pdbtools.make_pdblines_from_residues(residues,False)

    with open(args.output,'w') as outfile:
        outfile.write(''.join(headlines))
        outfile.write(''.join(reslines))
        outfile.write(''.join(tailLines))


main()
