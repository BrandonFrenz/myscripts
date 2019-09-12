#!/usr/bin/python
import argparse
import subprocess

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l','--list',help='The list of argo jobs to cancel')
    args = parser.parse_args()
    return args

def parse_list(listfile):
    jobs = []
    with open(listfile,'r') as lf:
        for line in lf:
            entry = line.split()
            #if entry[1] != 'Running':
            #    continue
            jobs.append(entry[0])
    return jobs

def cancel_jobs(jobs):
    for job in jobs:
        subprocess.call('argo delete {}'.format(job),shell=True)

def main():
    args = parseargs()
    jobs = parse_list(args.list)
    cancel_jobs(jobs)

main()

