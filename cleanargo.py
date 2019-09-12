#!/usr/bin/python
import subprocess
from collections import namedtuple
import re

Job = namedtuple('Job',['name','status','age'])

def load_list():
    jobs = subprocess.check_output(['argo','list'])
    jobdatum = jobs.splitlines()
    
    jobs = []
    for jobdata in jobdatum:
        data = jobdata.split()
        name = data[0]
        status = data[1]
        age = data[2]
        job = Job(name,status,age)
        jobs.append(job)
    return jobs

def clean_old_jobs(jobs):
    for job in jobs:
        if 'd' in job.age:
            age = re.sub('\D+','',job.age)
            if age > 1:
                if job.status != 'Running':
                    delete_job(job)

def delete_job(job):
    command = 'argo delete {}'.format(job.name)
    subprocess.call(command,shell=True)
    print command,job.name,job.status,job.age

def main():
    jobs =  load_list()
    clean_old_jobs(jobs)

main() 
