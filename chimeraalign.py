#!/usr/bin/env python
import argparse
import os, sys
import chimera
import Midas

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-ms', '--mobile_structure', help='The moving structure')
    parser.add_argument('-ss', '--static_structure', help='The structure staying in place')
    args = parser.parse_args()
    return args

def align_structures(static_structure, mobile_structure):
    chimera.openModels.open(static_structure)
    chimera.openModels.open(mobile_structure)
    opened = chimera.openModels.list()
    chimera.runCommand('mm #0 #1')
    Midas.write(opened[1], None, mobile_structure)
    Midas.close([opened[1]])


def main():
    args = parseargs()
    align_structures(args.static_structure, args.mobile_structure)

if __name__ == '__main__':
    main()
