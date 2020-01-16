#!/usr/bin/python
import argparse
import os
import subprocess

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--job_name')
    parser.add_argument('-b', '--bucket_name', default='argo-science-artifacts', help='The name of the bucket the file is in')
    args = parser.parse_args()
    return args

def download_files(jobname, bucket_name):
    jobid = '-'.join(jobname.split('-')[0:-1])
    cmnd = ['gsutil', 'cp', f'gs://{bucket_name}/{jobid}/{jobname}/*', './']
    subprocess.call(cmnd)
    print(cmnd)

def main():
    args = parseargs()
    download_files(args.job_name, args.bucket_name)

if __name__ == '__main__':
    main()
