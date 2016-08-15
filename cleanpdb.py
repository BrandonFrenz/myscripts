#!/usr/local/bin/python2.7
import sys

pdbs = []
count = 1
while count < len(sys.argv):
    pdbs.append(sys.argv[count])
    count +=1

for pdb in pdbs:
    newfile = []
    with open(pdb,'r') as pdbfile:
        for line in pdbfile:
            if line.startswith('ATOM'):
                restype = line[17:20]
                if 'UNK' in str(restype):
                    continue
                else:
                    newfile.append(line)

    newpdb = 'new'+str(pdb)
    with open(newpdb,'w') as newpdbfile:
        for line in newfile:
            newpdbfile.write(line)
