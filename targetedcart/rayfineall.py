#!/usr/bin/python
import argparse
import multiprocessing
import subprocess
import os


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure', help='The input structure')
    parser.add_argument('-d', '--density', help='The densit files')
    parser.add_argument('-i', '--iters', type=int, help='The number of iterations to run')
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
    command = ['sh', 'runlocal.sh', vars['structure'], '_{}'.format(vars['prefix']), vars['density']]
    subprocess.call(command)

def main():
    args = parseargs()
    pool = multiprocessing.Pool(args.cores)
    jobinputs = []
    for i in range(0, args.iters):
        jobinputs.append({'structure': args.structure, 'density': args.density, 'prefix': str(i+1)})
    print(jobinputs)
    pool.map(relax, jobinputs)

main()
