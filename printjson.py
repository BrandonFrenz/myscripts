#!/usr/bin/python
import argparse
import json

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--file',help='The json file to print')
    parser.add_argument('-o', '--output', help='The name of the output')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    with open(args.file,'r') as of:
        data = json.load(of)
    with open(args.output, 'w') as of:
        json.dump(data, of, indent=4)

if __name__ == '__main__':
    main()
