#!/usr/bin/env python
import argparse
import os
import subprocess

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-ms', '--mobile_structures', nargs="+", help='The the structures being aligned')
    parser.add_argument('-ss', '--static_structure', help='The structures everying is aligning to')
    args = parser.parse_args()
    return args

def launch_chimera(static_structure, mobile_structure):
    command = ['/Applications/Chimera.app/Contents/MacOS/chimera', '--nogui', '--nostatus', '--script', f'/Users/brandonfrenz/scripts/chimeraalign.py -ss {static_structure} -ms {mobile_structure}']
    subprocess.call(command)

def main():
    args = parseargs()
    for mobile_structure in args.mobile_structures:
        launch_chimera(args.static_structure, mobile_structure)

if __name__ == '__main__':
    main()
