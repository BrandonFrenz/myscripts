#!/usr/bin/python
import sys

beam = sys.argv[1]

bestrms = 99999
score = 0
with open(beam) as beamfile:
    for line in (beamfile):
        line = line.split()
        if len(line) == 3:
            rms = float(line[1])
            if rms < bestrms:
                bestrms = rms
                score = line[0]
print "best rms " +str(bestrms) + " score " + score
