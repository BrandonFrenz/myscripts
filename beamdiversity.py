#!/usr/bin/python
import argparse
import growerparser
import operator
import numpy as np

def main():
    args = parseargs()
    report_diversity(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--beamfile',help='The beamfile')
    parser.add_argument('-n','--beamnum',type=int,default=10,help='The number of beams to use for the averaging')
    args = parser.parse_args()
    return args

def report_diversity(args):
    beams = growerparser.parse_beamfile(open(args.beamfile,'r').readlines())
    beams = filter_beams(args,beams)
    print rms_variation(args,beams)

def filter_beams(args,beams):
    sorted_beams = sorted(beams, key=operator.attrgetter('score'))
    filtered_beams = []
    it = 0
    for beam in sorted_beams:
        filtered_beams.append(beam)
        if it >= args.beamnum:
            break
        it+=1
    return filtered_beams

def rms_variation(args,beams):
    caindex = 2
    ca_coms = []
    for beam in beams:
        for residue in beam.bbresidues:
            comx = 0
            comy = 0
            comz = 0
            for atom in residue:
                if atom.atomid == caindex:
                    comx+=atom.x
                    comy+=atom.y
                    comz+=atom.z
            ca_coms.append((comx,comy,comz))
    total_dist = 0
    total_cas = 0
    for beam in beams:
        resit = 0
        for residue in beam.bbresidues:
            for atom in residue:
                if atom.atomid == caindex:
                    com_coords = np.array(ca_coms[resit])
                    beam_coords = np.array((atom.x,atom.y,atom.z))
                    resit+=1
                    total_dist+= np.linalg.norm(com_coords-beam_coords)*2
                    total_cas+=1
    rms = np.sqrt(total_dist/total_cas)
    return rms


main()
