#!/usr/bin/python
import argparse
import os
import fileinput
import sys
import re

def main():
    args = parseargs()
    args = add_options(args)
    if args.mode == 'dump':
        dumpbeam_xml(args)
        dump_pdbs(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m','--mode',default='dump',help='The mode to run')
    parser.add_argument('-ex','--executable',default='~/Desktop/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease')
    parser.add_argument('-x','--xml',default='parallel.xml',help='The xml')
    parser.add_argument('-p','--pdb',default='input.pdb',help='the input pdb')
    parser.add_argument('-d','--density',default='density.mrc',help='The density map')
    parser.add_argument('-n','--native',default='',help='The native pdb')
    parser.add_argument('-b','--beams',default='na',help='The beam file to use')
    parser.add_argument('-tb','--taboobeams',default='na',help='The taboo beams to use')
    parser.add_argument('-s','--steps',default=0,help='The number of steps to run in the loop grower')
    args = parser.parse_args()
    return args

def add_options(args):
    args.options = '-default_max_cycles 200 -ignore_unrecognized_res -overwrite -mapreso 2 -beta -missing_density_to_jump'
    return args

def dump_pdbs(args):
    beams = 0
    taboo = 0
    if args.taboobeams != 'na':
        taboo = 1
    if args.beams != 'na':
        beams = 1
    native = ' '
    if args.native != '':
        native = ' -in::file::native %s '%(args.native)
    #exectuable, xml, density, native, beamread, beams, steps, tabooread, taboobeams
    command = ("%s -parser:protocol %s -s %s -edensity:mapfile %s%s-parser::script_vars readbeams=%s -parser::script_vars beams=%s -parser::script_vars steps=0 -parser::script_vars pcount=0 "
            "-parser::script_vars filterprevious=%s -parser::script_vars filterbeams=%s %s"""%(args.executable,args.xml,args.pdb,args.density,native,beams,args.beams,taboo,args.taboobeams,args.options))
    print command
    os.system(command)

def dumpbeam_xml(args):
    for line in fileinput.input(args.xml, inplace=1):
           sys.stdout.write(re.sub("dumpbeam=\d+", "dumpbeam=1", line))

main()
