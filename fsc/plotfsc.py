#!/usr/bin/env python
import argparse
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l1', '--log1', help='The log file with the fsc data in it')
    parser.add_argument('-l2', '--log2', help='The log file with the fsc data in it')
    parser.add_argument('-o', '--output', help='The output file name')
    args = parser.parse_args()
    return args

def parse_log(logfile):
    xs = []
    ys = []
    with open(logfile, 'r') as inf:
        for line in inf:
            data = line.split()
            x = float(data[1])
            y = float(data[-1])
            xs.append(x)
            ys.append(y)
    return xs, ys
 
def plot_xy(xs1, ys1, xs2, ys2, output):
    pp = PdfPages(output)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(r'FSC vs Resolution (1/$\AA$)')
    ax.plot(xs1, ys1, label='halfmap1')
    ax.plot(xs2, ys2, label='halfmap2')
    ax.set_xlabel(r'Resolution (1/$\AA$)')
    ax.set_ylabel('FSC')
    ax.legend(bbox_to_anchor=(1, 1), fancybox=True, shadow=True, ncol=5)
    plt.savefig(pp,format='pdf')
    pp.close()

def main():
    args = parseargs()
    xs1, ys1 = parse_log(args.log1)
    xs2, ys2 = parse_log(args.log2)
    plot_xy(xs1, ys1, xs2, ys2, args.output)

if __name__ == '__main__':
    main()
