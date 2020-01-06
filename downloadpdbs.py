#!/usr/bin/python
import argparse
import os

def main():
    args = parseargs()
    download_pdb(args.pdbids)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdbids',nargs="+",help='The pdb id of the structure to download')
    args = parser.parse_args()
    return args

def download_pdb(pdbs):
    for pdbid in pdbs:
        try:
            command = 'wget http://www.rcsb.org/pdb/files/%s.pdb.gz'%pdbid
            os.system(command) 
            os.system('gunzip %s.pdb.gz'%pdbid)
        except:
            print('error',pdbid,'could not be downloaded')

if __name__ == '__main__':
    main()
