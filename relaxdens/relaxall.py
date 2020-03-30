#!/usr/bin/python
import argparse
import multiprocessing
import subprocess
import os


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structures', nargs='+', help='The input structure files')
    parser.add_argument('-d', '--density', help='The densit files')
    parser.add_argument('-j', '--cores', type=int, default=1, help='The number of cores to use')
    args = parser.parse_args()
    return args

def relax(vars):
    rosetta = os.getenv('ROSETTA')
    if not os.path.isfile(vars['density']) or not os.path.isfile(vars['structure']):
        print('One or more input files missing')
        return
    if rosetta == None:
        print('ROSETTA environment variable not set')
        return
    command = ['{}/source/bin/rosetta_scripts.default.linuxgccrelease'.format(rosetta), '-s', vars['structure'], 
            '@relax_into_density.flags', '-parser:protocol', 'relax_into_density.xml', '-edensity:mapfile', vars['density'], '-database', '{}/database'.format(rosetta)]
    subprocess.call(command)
def main():
    args = parseargs()
    pool = multiprocessing.Pool(args.cores)
    jobinputs = []
    for structure in args.structures:
        jobinputs.append({'structure': structure, 'density': args.density})
    print(jobinputs)
    pool.map(relax, jobinputs)

main()
