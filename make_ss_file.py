#!/usr/bin/python
import argparse
import pdbtools

def main():
    args = parseargs()
    ssbonds = find_disulfs(args)
    write_ssbond_file(args,ssbonds)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',help='The pdbs')
    parser.add_argument('-d','--distance',type=float,default=5,help='The distance cutoff')
    parser.add_argument('-o','--output',help='The output disulf file')
    args = parser.parse_args()
    return args

def find_disulfs(args):
    residues = pdbtools.get_unopened_residue_list(args.pdb)
    connected = []
    count = 1
    ssbonds = []
    for res in residues:
        if res.name == 'CYS':
            for atom in res.atoms:
                if atom.atomid == ' SG ':
                    newcount = 1
                    for sres in residues:
                        if count == newcount:
                            continue
                        if sres.name == 'CYS':
                            for satom in sres.atoms:
                                if satom.atomid == ' SG ':
                                    distance = pdbtools.atom_dist(atom,satom)
                                    if distance < args.distance:
                                        ssbonds.append((count,newcount))
                                        if connected.count(count) > 2:
                                            print count
                                        if connected.count(newcount) > 2:
                                            print newcount
                                        if count in connected:
                                            continue
                                        if newcount in connected:
                                            continue
                                        connected.append(count)
                                        connected.append(newcount)
                        newcount+=1
        count+=1
    return ssbonds

def write_ssbond_file(args,ssbonds):
    with open(args.output,'w') as outfile:
        for ssbond in ssbonds:
            outfile.write(str(min(ssbond[0],ssbond[1]))+' '+str(max(ssbond[0],ssbond[1]))+'\n')

main()
