#!/usr/bin/env python
import argparse

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structures', nargs="+", help='The structures to analyze')
    parser.add_argument('-o', '--output', help='The output file')
    args = parser.parse_args()
    return args

def parse_disulf_results(structure_file):
    disulfides = []
    with open(structure_file, 'r') as inf:
        for line in inf:
            if line.startswith('SSBOND'):
                res1 = line.split()[3]
                res2 = line.split()[6]
                disulfides.append((res1,res2))
            if line.startswith('label'):
                disulf_index = line.split().index('dslf_fa13')
            if line.startswith('pose'):
                disulfide_score = float(line.split()[disulf_index])
    return disulfides, disulfide_score/len(disulfides)

def write_results(output, disulf_scores):
    with open(output, 'w') as of:
        of.write('Name, Disulfides, Average Disulfide Score\n')
        for struct in disulf_scores:
            of.write(struct)
            of.write(',')
            for dspair in disulf_scores[struct]['disulfides']:
                of.write(f'({dspair[0]} {dspair[1]}) ')
            of.write(',')
            of.write(str(disulf_scores[struct]['average disulfide score']))
            of.write('\n')

def main():
    args = parseargs()
    disulf_scores = {}
    for struct in args.structures:
        disulf_scores[struct] = {}
        disulfides, avg_score = parse_disulf_results(struct)
        disulf_scores[struct]['disulfides'] = disulfides
        disulf_scores[struct]['average disulfide score'] = avg_score
    write_results(args.output, disulf_scores)

if __name__ == '__main__':
    main()
