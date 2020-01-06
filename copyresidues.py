#!/usr/bin/python
import argparse
import re
import pdbtools

def main():
    args = parseargs()
    args = parse_residues(args)
    print len(args.res1),len(args.res2)
    copy_atoms(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p1','--pdb1',help='The first pdb')
    parser.add_argument('-p2','--pdb2',help='The second pdb')
    parser.add_argument('-r1','--res1',nargs="+",help='The first residues')
    parser.add_argument('-r2','--res2',nargs="+",help='The second residues')
    parser.add_argument('-o','--output',help='The name of the output model')
    args = parser.parse_args()
    return args

def parse_residues(args):
    residues = []
    for res in args.res1:
        res = res.split('-')
        chain = re.split('\d+',res[0])[0]
        lower = int(re.split('(\d+)',res[0])[1])
        upper = lower
        if len(res) > 1:
            upper = int(res[1])
        for i in range(lower,upper+1):
                residues.append((chain,i))
    args.res1 = residues
    
    residues = []
    for res in args.res2:
        res = res.split('-')
        chain = re.split('\d+',res[0])[0]
        lower = int(re.split('(\d+)',res[0])[1])
        upper = lower
        if len(res) > 1:
            upper = int(res[1])
        for i in range(lower,upper+1):
                residues.append((chain,i))
    args.res2 = residues
    
    return args

def copy_atoms(args):
    residues1 = pdbtools.get_unopened_residue_list(args.pdb1)
    residues2 = pdbtools.get_unopened_residue_list(args.pdb2)

    resmatch = 0
    resis_to_copy = []
    for residue in residues1:
        print (residue.chain,residue.num)
        if (residue.chain,residue.num) in args.res1:
            resis_to_copy.append(residue)

    copiedres = 0
    updated_residues = []
    for residue in residues2:
        if (residue.chain,residue.num) in args.res2:
            atoms = []
            for it in range(0,len(residue.atoms)):
                atom1 = resis_to_copy[copiedres].atoms[it]
                atom2 = residue.atoms[it]
                atom2.x = atom1.x
                atom2.y = atom1.y
                atom2.z = atom1.z
                atoms.append(atom2)
            residue.atoms = atoms
            copiedres+=1
        updated_residues.append(residue)
    pdbtools.write_resis_to_pdb(updated_residues,args.output)


main()
