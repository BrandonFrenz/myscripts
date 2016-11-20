#!/bin/usr/python2.7/
import glob
import os
import sys

pdblist = []
i=1
while i < len(sys.argv):
    pdblist.append(sys.argv[i])
    i+=1
for pdb in pdblist:
    print '<Template pdb="%s" weight=1.0 cst_file="AUTO"/>'%pdb
