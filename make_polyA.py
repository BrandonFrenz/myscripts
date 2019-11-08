#!/usr/bin/env python
import argparse
from pyrosetta import *
from pyrosetta.rosetta import core
init()

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure', help='Pose from file')
    parser.add_argument('-o', '--output', default='polya.pdb', help='The output file name')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    pose = pose_from_file(args.structure)
    res = 1
    while res <= pose.size():
        if pose.residue(res).type().is_disulfide_bonded():
            ds_partner = core.conformation.get_disulf_partner(pose.conformation(), res)
            core.conformation.break_disulfide(pose.conformation(), res, ds_partner)
        pyrosetta.toolbox.mutants.mutate_residue(pose, res, 'A')
        res+=1
    pose.dump_pdb(args.output)
    print(pose.size())

if __name__ == '__main__':
    main()
