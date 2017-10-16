#!/usr/bin/python
import subprocess
import sys
import os

def main():
    pdbs = []
    it = 1
    while it < len(sys.argv):
        pdbs.append(sys.argv[it])
        it+=1
    command = 'pymol %s &>/dev/null'%' '.join(pdbs)
    FNULL = open(os.devnull, 'w')
    subprocess.call(command,shell=True,stdout=FNULL)

main()
