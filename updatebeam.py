#!/usr/bin/python2.7
import sys

beams = []
count = 1
while count < len(sys.argv):
    beams.append(sys.argv[count])
    count+=1
for beam in beams:
    with open(beam,'r') as oldbeam:
        newbeam = []
        for line in oldbeam:
            line = line.split()
            if len(line) == 2:
                line.append("0.5")
            line+="\n"
            newline = " ".join(line)
            newbeam.append(newline)
        with open(beam,'w') as newbeamfile:
            for  i in newbeam:
                newbeamfile.write(i)
