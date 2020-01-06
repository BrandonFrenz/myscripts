#!/usr/bin/python
import subprocess
import sys
import os

def main():
    args = []
    it = 1
    while it < len(sys.argv):
        args.append(sys.argv[it])
        it+=1
    command = '/Applications/Pymol.app/Contents/bin/pymol %s &'%' '.join(args)
    FNULL = open(os.devnull, 'w')
    subprocess.call(command,shell=True,stdout=FNULL)

main()
