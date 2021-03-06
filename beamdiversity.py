#!/usr/bin/python
import argparse
import growerparser
import operator
import numpy as np

def main():
    args = parseargs()
    print get_diversity(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--beamfile',help='The beamfile')
    parser.add_argument('-n','--beamnum',type=int,default=10,help='The number of beams to use for the averaging')
    args = parser.parse_args()
    return args

def get_diversity(args):
    beams = []
    with open(args.beamfile,'r') as beamfile:
        beams = growerparser.parse_beamfile(beamfile.readlines())
    beams = filter_beams(args,beams)
    rms = rms_variation(args,beams)
    return rms

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
    totalres = len(beams[0].bbresidues)
    ca_coms = [(0,0,0)]*totalres
    totalbeams = len(beams)
    for beam in beams:
        resit = 0
        for residue in beam.bbresidues:
            comx = ca_coms[resit][0]
            comy = ca_coms[resit][1]
            comz = ca_coms[resit][2]
            for atom in residue:
                if atom.atomid == caindex:
                    comx+=atom.x/totalbeams
                    comy+=atom.y/totalbeams
                    comz+=atom.z/totalbeams
            ca_coms[resit] = ((comx,comy,comz))
            #if resit == 166:
            #    print comx,comy,comz
            #    exit()
            resit+=1
    total_dist = 0
    total_cas = 0
    for beam in beams:
        resit = 0
        for residue in beam.bbresidues:
            for atom in residue:
                if atom.atomid == caindex:
                    com_coords = np.array(ca_coms[resit])
                    beam_coords = np.array((atom.x,atom.y,atom.z))
                    dist = np.linalg.norm(com_coords-beam_coords)*2
                    if dist > 0:
                        total_cas+=1
                        resit+=1
                        total_dist+=dist
    rms = np.sqrt(total_dist/total_cas)
    return rms


#main()
