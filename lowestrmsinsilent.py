#!/bin/python2.7/
import sys

def main():

    #load silentfiles
    silentfiles = []
    it=1
    while it < len(sys.argv):
        silentfiles.append(sys.argv[it])
        it+=1
    
    #read and get all the structures below rmscutoff
    rmscutoff = 5.0
    count = 0
    for sf in silentfiles:
        count+=1
        newfile = []
        with open(sf, 'r') as silentfile:
            header = True
            firstscore = True
            rms = 999
            for line in silentfile:
                if header == True:
                    newfile.append(line)
                    if "REMARK" in line:
                        header = False
                    continue
                if header == False and line.startswith("SCORE:"):
                    splitline = line.split()
                    rms = float(splitline[5])
                if rms < rmscutoff:
                    newfile.append(line)
        writenewfile(newfile,count)

def writenewfile(filelist,count):
    if len(filelist) < 4:
        return
    with open("nn%s.silent"%count,"ab+") as newfile:
        for line in filelist:
            newfile.write(line)

main()
