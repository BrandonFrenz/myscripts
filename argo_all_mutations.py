#!/usr/bin/python
import argparse
import subprocess
import re
from collections import defaultdict
import os
import json
from Bio.PDB.PDBParser import PDBParser
import amino_acids


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--structure_file',help='The pdbfile')
    parser.add_argument('-bs','--benchmark_set',help='The benchmark set')
    parser.add_argument('-b','--bucket_name',default='argo-science-artifacts',help='The bucket name')
    parser.add_argument('-mf','--mutfile_key',help='The location to store all the mutation.json files')
    parser.add_argument('-w','--workflow',default='~/argo-workflows/workflows/cartesian-ddg/workflow.yaml',help='The argo workflow')
    parser.add_argument('-o','--output',default='brandon/',help='The output prefix on gce')
    parser.add_argument('-p','--pdbkey',help='The pdb key to use. This script assumes the pdb has been relaxed and is on gcs')
    parser.add_argument('-m','--missing',default=None,help='The file which contains missing wild types number pairs')
    parser.add_argument('-u','--upload',action='store_true',default=False,help='Upload the mutations files or not?')
    args = parser.parse_args()
    return args

def pdb_to_all_mutations_file(pdbfile,missing):
    wildtypes = []
    structure = parse_pdb(pdbfile)
    model = structure[0]
    for chain in model:
        for res in chain:
            wildtypes.append(amino_acids.longer_names[res.get_resname()])

    if not os.path.isdir('mutfiles/'):
        os.mkdir('mutfiles/')
    mutfiles = []
    for num,wt, in enumerate(wildtypes):
        if missing != None and [wt,num] not in missing:
            continue
        mutstrings = []
        for mut in amino_acids.one_letter_names:
            mutstring = ','.join([wt,str(num+1),mut])
            mutstrings.append(mutstring)
        mutfilename = 'mutfiles/mutfile_{}_{}.json'.format(num,wt)
        with open(mutfilename,'w') as mf:
            mf.write(json.dumps(mutstrings,indent=4))
        mutfiles.append(mutfilename)
    return mutfiles


def parse_pdb(pdbfile):
    parser = PDBParser()
    base = os.path.basename(pdbfile)
    structure = parser.get_structure(base,pdbfile)
    return structure


def upload_mutations(filenames,bucket_name,bucket_key):
    #basenames = []
    #for filename in filenames:
    #    basename = os.path.basename(filename)
    #    basenames.append(basename)
    #if bucket_key.endswith('/'):
    #    bucket_key = bucket_key.strip('/')
    sutilcmnd = 'gsutil -m cp {} gs://{}/{}/'.format(" ".join(filenames),bucket_name,bucket_key)
    print sutilcmnd
    subprocess.call(sutilcmnd,shell=True)

def run_argo(mutfiles,workflow,bucket_name,pdbkey,mutation_prefix,output_prefix,relax_input):
    for mutfile in mutfiles:
        mutbase = os.path.basename(mutfile)
        mutkey = os.path.join(mutation_prefix,mutbase)
        pdb_tag = re.sub('.json','',mutbase)
        outfile = pdb_tag+'_results.json'
        output_key = os.path.join(output_prefix,outfile)
        command = ('argo submit {} -p input-model="{}" -p mutation-file="{}" -p relax-input="{}" -p bucket-name="{}" -p output-bucket="{}" -p output-key="{}" -p use-custom-output="true"'
        .format(workflow,pdbkey,mutkey,relax_input,bucket_name,bucket_name,output_key))
        subprocess.call(command,shell=True)

def parse_missing(missing_file):
    if missing_file == None:
        return None
    with open(missing_file,'r') as mf:
        data = json.load(mf)
    return data

def main():
    args = parseargs()
    missing = parse_missing(args.missing)
    filenames = pdb_to_all_mutations_file(args.structure_file,missing)
    print len(filenames)
    relax_input = 'false'
    if args.upload:
        upload_mutations(filenames,args.bucket_name,args.mutfile_key)
    run_argo(filenames,args.workflow,args.bucket_name,args.pdbkey,args.mutfile_key,args.output,relax_input)

main()
