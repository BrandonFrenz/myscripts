#!/usr/bin/python
import argparse
import growerparser
import operator

def main():
    args = parseargs()
    res_scores = get_average_res_score(args)
    print_res_scores(args,res_scores)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--beams',nargs="+",help='The beam files')
    parser.add_argument('-o','--outfile',default='average_res_score.txt',help='The name of the output file')
    args = parser.parse_args()
    return args

def get_average_res_score(args):
    res_scores = {}
    for beam in sorted(args.beams):
       beams = growerparser.parse_unopened_beamfile(beam)
       nres = beams[0].total_res
       best_score = sorted(beams, key=operator.attrgetter('score'))[0].score
       average_res = best_score/nres
       res_scores[nres] = average_res
    return res_scores

def print_res_scores(args,res_scores):
    with open(args.outfile,'w') as outfile:
        for key in sorted(res_scores.keys()):
            outfile.write(str(key)+','+str(res_scores[key])+'\n')

main()
