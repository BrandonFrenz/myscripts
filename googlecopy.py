#!/usr/bin/python
import argparse
import os
import subprocess

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--job_name')
    args = parser.parse_args()
    return args

def download_files(jobname):
    jobid = '-'.join(jobname.split('-')[0:-1])
    cmnd = 'gsutil cp gs://argo-science-artifacts/{}/{}/* ./'.format(jobid,jobname)
    subprocess.call(cmnd,shell=True)
    print(cmnd)

def main():
    args = parseargs()
    download_files(args.job_name)

if __name__ == '__main__':
    main()
