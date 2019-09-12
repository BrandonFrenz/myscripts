#!/usr/bin/env python
import argparse

def parseargs():
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    return args

def main():
    args = parseargs()

if __name__ == '__main__':
    main()
