#!/usr/bin/python
import argparse
import os
import multiprocessing

def main():
    args = parseargs()

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--cores',help='The number of cores to use')
    args = parser.parse_args()
    return args

def run(args):
    command = 'sh run.sh %s'%args
    os.system(command)

def run_parallel(args)
    jobcount = 0
    totaljobs = len(args.inputs)
    while jobcount < totaljobs:
        if __name__ == '__main__':
            jobs = []
            for i in range(0,args.cores):
                jobinut = args.inputs[jobcount]
                p = multiprocessing.Process(target=extractfile, args=(jobinput),)
                jobs.append(p)
                p.start() 
                jobcount = jobcount + 1
                if jobcount >= totaljobs:
                    break
            for i in jobs:
                i.join()

main()
