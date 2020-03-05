#!/usr/bin/python
import argparse
import multiprocessing
import subprocess

EXECUTABLE = '/home/bfrenz/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease'
DATABASE = '/home/bfrenz/Rosetta/main/database'
def relax(structure):
    command = [EXECUTABLE, '-s', structure, 
            '@relax_into_density.flags', '-parser:protocol', 'relax_into_density.xml', '-edensity:mapfile', 'nosymm.mrc', '-database', DATABASE]
    subprocess.call(command)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structures', nargs='+')
    parser.add_argument('-j', '--cores', type=int, default=1, help='The number of cores to use')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    pool = multiprocessing.Pool(args.cores)
    pool.map(relax, args.structures)

main()
