#!/usr/bin/python
#original script written by Dan Farrel modified by Brandon Frenz
import fileinput
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import argparse
import re

def main():
    args = parseargs()
    #if len(args.silentfiles) > 0:
    #    sfilerms = readsilentfilerms(args.silentfiles,'sfilerms.txt')
    #    makeplot(sfilerms,'plot')
    #else:
    plotlist(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--silentfiles',nargs="+",help='The silentfiles')
    parser.add_argument('-l','--log',help='The log file')
    parser.add_argument('-r','--roundcount',help='Which round results do you want?')
    parser.add_argument('-o','--out',help='Where to put the outputfiles')
    parser.add_argument('-sd','--sfiledata',nargs="+",help='the silentrms data')
    args = parser.parse_args()
    return args

# begin opening silent files for plotting silent rms
def readsilentfilerms(silentfiles,name):
    silentrmsdata = []
    for x in silentfiles:
        # don't want any assembled files so we if continue them out
        if x.startswith("assembled"):
            continue
        for line in fileinput.input(x):
            # now we get the data but first get rid of the header file (score score)
            if line.startswith("SCORE:"):
                if line.startswith("SCORE:     score"):
                    continue
                res = line.rsplit("_")[-2]
                rms = float(line.rsplit()[-2])
                if rms == 0:
                    return
                silentrmsdata.append([rms, res])
    # this is code to sort the data, not super useful, but it's there if you want it
    silentrmsdata.sort()
    strsilentrms = []
    with open('%s/%s'%(args.out,name),'w') as rmsdatafile:
        for x in silentrmsdata:
            # write with 3 before . and 4 after . I don't get what the 0 and 1 mean, but I think it has to do with string stuff
            strsilentrms.append("%3.4f,%s\n" % (x[0], x[1]))
            rmsdatafile.write("%3.4f,%s\n" % (x[0], x[1]))
    return strsilentrms

# beginning of matplotlib graphing/printing first is a combined graph
# get our data from the outputted text files as from above (may be useful in the future to avoid this
# by just getting the data straight from the lists above, but this is easy because we almost always want to
# move the txt files by hand for various comparisons so this works out just fine
def makeplot(sfilename,name):
    silentrmsplot = np.genfromtxt(sfilename, delimiter=',', dtype=float, names=['yd', 'xd'])
    
    # ###this one will print silent rms plot######
    # begin with figure
    plt.figure()
    # set our axis
    ax2 = plt.gca()
    # plot from our silentrmsplot
    ax2.plot(silentrmsplot['xd'], silentrmsplot['yd'], 'o', c='r', alpha=1, markeredgecolor='none')
    # set yscale to log (again, can be missleading
    ax2.set_yscale('log')
    # set our yticks and add a formatter because matplotlib is weird
    ax2.set_yticks([0, 1, 2, 5, 10, 20, 50])
    for axis in [ax2.xaxis, ax2.yaxis]:
        axis.set_major_formatter(ScalarFormatter())
    # set our ylimits because matplotlib really wants to show data below 0??
    ax2.set_ylim(0, max(silentrmsplot['yd'])*1.1)
    # this is in order to set the limits in the x axis in the plot
    # this is only useful in a log plot because the data bunched up at the top looks weird if it gets
    # cut off halfway
    ax2.set_xlim([min(silentrmsplot['xd']) - 5, max(silentrmsplot['xd']) + 5])
    # this will add the line accross the screen at 2
    ax2.axhline(y=2, linewidth=1, linestyle='--')
    # save our figure
    plt.savefig('%s'%name, bbox_inches='tight')
   
def plotlist(args):
    for data in args.sfiledata:
        name = re.split('.txt',data)[0]
        name += ".png"
        makeplot(data,name)

main()
