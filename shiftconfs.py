#!/usr/bin/env python
import argparse

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--confs', help='The conformers file')
    parser.add_argument('-o', '--output', help='The name of the output file')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    with open(args.confs, 'r') as inf:
        lines = inf.readlines()
    newlines = []
    for line in lines:
        if line.startswith('HETATM'):
            atomid = ''.join(line[12:16])
            newatomid = ' '+''.join(atomid[0:-1])
            newline = ''.join(line[0:12])+newatomid+''.join(line[16:])
            print(line)
            print(newline)
            #print("'{}'".format(atomid))
            #print("'{}'".format(newatomid))S
            newlines.append(newline)
        else:
            newlines.append(line)
    with open(args.output, 'w') as of:
        of.write(''.join(newlines))

if __name__ == '__main__':
    main()
