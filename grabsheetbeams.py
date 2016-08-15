#!/usr/bin/python
import argparse

def main():
    args = parseargs()
    get_sheet_beams(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--beamfiles',nargs="+",help='The beamfiles to check')
    args = parser.parse_args()
    return args

def get_sheet_beams(args):
    for beamfile in args.beamfiles:
        sheetbeams = []
        bf = open(beamfile,'r').readlines()
        first_line = bf[0]
        currentbeam = []
        for line in bf:
            if line == first_line:
                if len(currentbeam) >5:
                    sheetbeams.append(currentbeam)
                currentbeam = []
                currentbeam.append(line)
            else:
                currentbeam.append(line)
        if len(currentbeam) > 5:
            sheetbeams.append(currentbeam)
        filename = 'sheet'+beamfile
        with open(filename,'w') as sbfile:
            for sb in sheetbeams:
                for sbline in sb:
                    sbfile.write(sbline)
                



main()
