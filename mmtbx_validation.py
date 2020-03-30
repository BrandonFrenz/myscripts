#!/usr/bin/env python
import argparse

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structure', help='The structure to analyze')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()

if __name__ == '__main__':
    main()
