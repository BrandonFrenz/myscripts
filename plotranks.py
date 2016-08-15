#!/usr/local/bin/python2.7
import numpy as np
import scipy
import sys
from matplotlib import pyplot as plt
import glob
from matplotlib.backends.backend_pdf import PdfPages


with PdfPages('multipage_pdf.pdf') as pdf:
    with open('structureids.txt', "r") as ids:
        for structureid in ids:
            fig = plt.figure()
            ax = plt.subplot(111)
            cleanid = structureid.replace('\n','')
            plt.title('%s'  % cleanid)
            hasdata = False
            for experiment in glob.glob("*%s*" %cleanid):
                hasdata = True
                ranks = open('%s' % experiment, "r")
                cycles = []
                rankings = []
                mostresidues = 0
                cyclecount = 0
                rms = 0
                for line in ranks:
                    split = line.split()
                    cycle = split[1]
                    if '.0' in cycle:
                        if float(cycle) > cyclecount:
                            cyclecount = float(cycle)
                            rms = split[8]
                            rank = split[10]
                            if float(rms) > 1.5:
                                rank = 300
                            cycles.append(cycle)
                            rankings.append(rank)
                if len(cycles) > mostresidues:
                    mostresidues  = len(cycles)
                ax.axis([0, mostresidues, 0, 300])
                legendname = "%s" % experiment
                legendname = legendname.replace('%s.txt' % cleanid, "")
                color = 'r-'
                if 'sheetsampler' in legendname:
                    color = 'b-'
                legendname = legendname + "_" + str(rms)
                ax.plot(cycles, rankings, color, label=legendname )
            if hasdata:
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height * 0.1,
                    box.width, box.height * 0.9])

                ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                    fancybox=True, shadow=True, ncol=5)
                pdf.savefig()
                #plt.show()
                plt.close()
print 'results in multipage_pdf.pdf'
