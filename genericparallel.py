#!/usr/bin/python
import argparse
import multiprocessing
import subprocess

#given a set of pdbs a command line argument and a number of cores this script will run jobs for each pdb in parallel

def main():
    args = parseargs()
    run_parallel(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--command',help='The command to run')
    parser.add_argument('-p','--pdbs',nargs="+",help='The pdbs to run on')
    parser.add_argument('-j','--cores',type=int,default=1,help='The number of computers to use')
    args = parser.parse_args()
    return args

def run_parallel(args):
    if __name__ == '__main__':
        jobs = []
        for pdb in args.pdbs:
            jobs.append((args,pdb))
        p = multiprocessing.Pool(args.cores,maxtasksperchild=1)
        p.map(run, jobs)

def run(inputs):
    args = inputs[0]
    pdb = inputs[1]
    command = args.command+' '+pdb
    subprocess.call(command,shell=True)

    

main()
