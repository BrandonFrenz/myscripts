#!/usr/bin/python
import sys

def main():
    beams = []
    it = 1;
    while it < len(sys.argv):
        beams.append(sys.argv[it])
        it+=1
    
    for beamI in beams:
        with open(beamI) as beamfile:
            beam = []
            beamlist = []
            for line in beamfile:
                line = line.split()
                if len(line) == 5:
                    if len(beam) != 0:
                        beamlist.append(beam)
                        beam = []
                beam.append(line)
            beamlist.append(beam)
        beamlist = removeredundant(beamlist)
        writebeamlist(beamlist, beamI)

def writebeamlist(beamlist, filename):
    with open(filename, 'ab+') as newbeamfile:
        it = 0
        for i in beamlist:
            for ii in beamlist[it]:
                newbeamfile.write(" ".join(ii)+"\n")
            it+=1

def removeredundant(beamlist):
    newbeamlist = []
    for beam in beamlist:
        addbeam = True
        for newbeam in newbeamlist:
            if beam[0] == newbeam[0]:
                if beam[1] == newbeam[1]:
                    addbeam = False
        if addbeam == True:
            newbeamlist.append(beam)
    return newbeamlist

main()
