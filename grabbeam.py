#!/usr/local/bin/python2.7
import argparse

def main():
    args = parseargs()
    currentbeam = findmatch(args)
    printbeam(currentbeam)
    

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rms", type=float, help="the rms of the structure you are looking for")
    parser.add_argument("-b", "--beams", help="the beamfile to search")
    args = parser.parse_args()
    return args

def findmatch(args):
    with open(args.beams,"r") as beams:
        currentbeam = []
        matchingbeam = []
        mindiff = -1
        for line in beams:
            if len(line.split()) == 5:
                currentbeam = []
            currentbeam.append(line)
            if len(line.split()) == 3:
                rms = float(line.split()[1])
                diff = abs(rms-args.rms)
                if diff < mindiff or mindiff == -1:
                    mindiff = diff
                    matchingbeam = currentbeam

        return matchingbeam

def printbeam(currentbeam):
    for i in currentbeam:
        print i.strip()

main()
