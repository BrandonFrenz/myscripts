#!/usr/bin/python
import argparse
import amino_acids
import pdbtools

def main():
    args = parseargs()
    move_to_end(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',help='The pdb file')
    parser.add_argument('-o','--output',help='The name of the output')
    args = parser.parse_args()
    return args

def move_to_end(args):
    #Grab the head and tail lines before we manipulate the residues
    header = True
    tailer = True
    headlines = []
    taillines = []
    with open(args.pdb,'r') as pdbfile:
        pdblines = pdbfile.readlines()
        it = 0
        while header == True:
            if not (pdblines[it].startswith('ATOM') or pdblines[it].startswith('HETATM')):
                headlines.append(pdblines[it])
            else:
                header = False
            it+=1
        it = len(pdblines)-1
        while tailer == True:
            if not (pdblines[it].startswith('ATOM') or pdblines[it].startswith('HETATM')):
                taillines.append(pdblines[it])
            else:
                tailer = False
            it-=1

    residues = pdbtools.get_unopened_residue_list(args.pdb)
    movedresidues = []
    for res in residues:
        if res.name not in amino_acids.longer_names:
            continue
        movedresidues.append(res)
    for res in residues:
        if res.name not in amino_acids.longer_names:
            movedresidues.append(res)
    reslines = pdbtools.make_pdblines_from_residues(movedresidues,False)
    with open(args.output,'w') as outfile:
        outfile.write(''.join(headlines))
        outfile.write(''.join(reslines))
        outfile.write(''.join(taillines))

main()
