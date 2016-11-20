#!/usr/bin/python
import argparse

def main():
    args = parseargs()
    make_csv(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',help='the input file')
    parser.add_argument('-o','--output',help='the outputfile')
    args = parser.parse_args()
    return args

def make_csv(args):
    file_lines = []
    with open(args.input,'r') as infile:
        file_lines = infile.readlines()
    with open(args.output,'w') as outfile:
        for line in file_lines:
            elements = line.split()
            outfile.write(','.join(elements)+"\n")

main()
