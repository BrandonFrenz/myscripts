#!/usr/bin/python
import argparse
import glob
import numpy

def main():
    args = parseargs()
    ddgdata = grabscores(args)
    ddgstats = calculatestatistics(args,ddgdata)
    reportstats(ddgstats)
    #writeddg(ddgdata)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdbs',nargs="+",help='the pdbs to match')
    args = parser.parse_args()
    return args

def grabscores(args):
    allddgdata = {}
    for pdb in args.pdbs:
        pdbid = pdb.split('.pdb')[0]
        scfiles = glob.glob('%s*.sc'%pdbid)
        perpdbddg = {}
        for scfile in scfiles:
            perresddg = {}
            with open(scfile,'r') as scfile:
                pres = 1
                first = True
                for line in scfile:
                    data = line.split()
                    res = int(data[0])
                    subres = data[2]
                    ddg = data[3]
                    if res != pres and first == False:
                        perpdbddg[pres] = perresddg
                        perresddg = {}
                    first = False
                    perresddg[subres] = ddg
                    pres = res
                perpdbddg[pres] = perresddg
        allddgdata[pdbid] = (perpdbddg)
    return allddgdata

def calculatestatistics(args,ddgdata):
    datasummary = {}
    pdbid = args.pdbs[0].split('.pdb')[0]
    totalres = len(ddgdata[pdbid])
    it=1
    firstpdbid = args.pdbs[0].split('.pdb')[0]
    residues = ddgdata[firstpdbid].keys()
    firstres = residues[0]
    mutations = ddgdata[firstpdbid][firstres].keys()
    ddgstats = {}
    for res in residues:
        mutdic = {}
        for mut in mutations:
            mutddgs = []
            for pdb in ddgdata.keys():
                mutddgs.append(float(ddgdata[pdb][res][mut]))
            stddev = numpy.std(mutddgs)
            mean = numpy.mean(mutddgs)
            stats = []
            stats.append(mean)
            stats.append(stddev)
            mutdic[mut] = stats
        ddgstats[res] = mutdic
    return ddgstats

def reportstats(ddgstats):
    with open('ddgstats.csv','w') as statsfile:
        residues = ddgstats.keys()
        mutations = ddgstats[residues[0]].keys()
        header = 'residue,'
        for mut in mutations:
            header += '%s mean,%s stddev,'%(mut,mut)
        header += "\n"
        statsfile.write(header)
        for res in ddgstats.keys():
            reportline = str(res)
            for mut in ddgstats[res].keys():
                mean = ddgstats[res][mut][0]
                stddev = ddgstats[res][mut][1]
                mean = "{0:.4f}".format(mean)
                stddev = "{0:.4f}".format(stddev)
                reportline += ','+str(mean)+","+str(stddev)
            reportline +='\n'
            statsfile.write(reportline)

def writeddg(args,ddgdata):
    print "writing ddg"
    for pdb in ddgdata.keys():
        with open('%s_ddg.csv'%pdb,'w') as ddgdatafile:
            header = "res,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y\n"
            ddgdatafile.write(header)
            it = 1
            for res in sorted(ddgdata[pdb].keys()):
                resdata = str(res)
                for mut in sorted(ddgdata[pdb][res].keys()):
                    ddg = ddgdata[pdb][res][mut]
                    resdata+=",%s"%str(ddg)
                ddgdatafile.write(resdata+"\n")


main()
