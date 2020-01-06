#!/usr/bin/python
import argparse
import re
import json
from collections import defaultdict
from Bio.PDB import PDBList
import os

def main():
    args = parseargs()
    entries = parse_hhr(args)
    isolate_alignments(entries)
    if args.dssp:
        turn_off_nodssp(entries)
    if args.mode == 'get pdbs':
        get_pdbs(args,entries)
    if args.mode == 'top N':
        get_topN(args,entries)
        write_entries(args,entries)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-hhr','--hhr',help='The hhr file')
    parser.add_argument('-m','--mode',help='What to do')
    parser.add_argument('-n','--topN',type=int,default=0,help='How many models to download')
    parser.add_argument('-o','--output',default='out.hhr',help='The output hhr')
    parser.add_argument('-af','--alignment_format',default='grishin',help='The format of the output alignment')
    parser.add_argument('-dssp','--dssp',default=True,help='Turn off print the alignments with no dssp')
    args = parser.parse_args()
    return args

def parse_hhr(args):
    entries = defaultdict(dict)
    with open(args.hhr,'r') as hhrfile:
        hitbox = False
        first = True
        enum = 0
        elines = []
        header = True
        headerlines = []
        hasdssp = False
        ranklines = []
        ranks = False
        for line in hhrfile:
            if line.startswith('Query'):
                target = line.split()[1]
                entries['target'] = target.strip(',')
            if hitbox == True and not line.strip() == '' and 'No' not in line:
                num = int(line.split()[0])
                pdbid = line.split()[1]
                entries[num]['pdbid'] = pdbid
            if line.strip().startswith('No Hit'):
                ranks = True
                hitbox = True
                header = False
                entries['header']['lines'] = headerlines
            if header:
                headerlines.append(line)
            if line.startswith('No') and 'Hit' not in line and len(line.split()) == 2:
                hitbox = False
                if not first:
                    entries[enum]['lines'] = elines
                    entries[enum]['dssp'] = hasdssp
                    entries[enum]['print'] = True
                    elines = []
                    hasdssp = False
                entries['ranks'] = ranklines
                ranks = False
                first = False
            if ranks == True:
                ranklines.append(line)
            if not first:
                if line.startswith('No'):
                    enum = int(line.split()[1])
                elines.append(line)
            if 'dssp' in line:
                hasdssp = True
        entries[enum]['lines'] = elines
        entries[enum]['print'] = True
        entries[enum]['dssp'] = hasdssp
    return entries 

#prints the entries to the screen
def print_entries(args,entries):
    print json.dumps(entries,sort_keys=True,indent=4,separators=(',',': '))

#writes the entires to the outputfile in proper .hhr format.
def write_entries(args,entries):
    with open(args.output,'w') as outfile:
        if args.alignment_format == 'hhr':
            outfile.write(''.join(entries['header']['lines']))
            outfile.write(''.join(entries['ranks']))
            for entry in entries:
                try:
                    if entries[entry]['print']:
                        outfile.write(''.join(entries[entry]['lines']))
                except:
                    continue
            outfile.write('Done!')
        elif args.alignment_format == 'grishin':
            for entry in entries:
                if entries[entry]['print'] != True:
                    continue
                outfile.write('## 1XXX_ '+entries[entry]['pdbid']+'\n')
                outfile.write('# hhsearch\n')
                outfile.write('scores_from_program: 0 1.00\n')

def get_pdbs(args,entries):
    for entry in entries.keys():
        if entry == 'header':
            continue
        if args.topN != 0 and entry > args.topN:
            continue
        try:
            pdbid = entries[entry]['pdbid'].split('_')[0]
            download_pdb(pdbid)
        except:
            print 'No pdb for',entry

def download_pdb(pdbid):
    pdbl = PDBList()
    pdbfile = "pdb"+pdbid+".ent"
    if not os.path.isfile('pdbs/%s'%pdbfile):
        try:
            pdbl.retrieve_pdb_file(pdbid,pdir='pdbs')
        except:
            print 'could not download', pdbid

def isolate_alignments(entries):
    for entry in entries:
        try:
            continue
            isint = int(entry)
            alignment = []
            for line in entries[entry]['lines']:
                try:
                    if entries['target'] in line.split()[1]:
                        alignment.append(line.split()[3])
                except:
                    continue
            entry['alignment'] = ''.join(alignment)
        except:
            continue


#accepts only models who's rank is greater than the specified value
def get_topN(args,entries):
    for entry in entries.keys():
        try:
            intentry = int(entry)
            if intentry > args.topN:
                entries[entry]['print'] = False
        except:
            continue

#Turns off the printing of any alignment which doesn't have dssp
def turn_off_nodssp(entries):
    for entry in entries:
        try:
            if entries[entry]['dssp'] == False:
                entries[entry]['print'] = False
        except:
            continue

main()
