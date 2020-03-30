#!/usr/bin/env python
import argparse
import mmtbx
from mmtbx.validation.molprobity import molprobity
#from validation import molprobity
import iotbx.pdb.hierarchy
from mmtbx.command_line import load_model_and_data
from mmtbx import model

master_phil_str = """
clashscore = True
  .type = bool
ramalyze = True
  .type = bool
omegalyze = True
  .type = bool
rotalyze = True
  .type = bool
cbetadev = True
  .type = bool
nqh = True
  .type = bool
rna = True
  .type = bool
model_stats = True
  .type = bool
restraints = True
  .type = bool
rfactors = True
  .type = bool
real_space = True
  .type = bool
waters = True
  .type = bool
seq = True
  .type = bool
xtriage = False
  .type = bool
"""

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure', help='The input structure')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()

    pdb_in = iotbx.pdb.hierarchy.input(file_name=args.structure)
    model = mmtbx.model.manager(pdb_in.input)
    #print(pdb_in.hierarchy.only_model())
    #print(pdb_in.hierarchy.atom_size())
    #pdb_in.construct_hierarchy()
    #pdb_hierarchy = mmtbx.secondary_structure.get_pdb_hierarchy(args.structure)
    #print(type(pdb_in.hierarchy))
    print(molprobity(model))

if __name__ == '__main__':
    main()
