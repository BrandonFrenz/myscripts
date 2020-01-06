#!/usr/bin/python
import argparse
from collections import defaultdict
from Bio.PDB import PDBList
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import os
import multiprocessing
import random
import amino_acids
from copy import deepcopy
import pdbtools
import glob
import re
import time
import matplotlib
import matplotlib.pyplot as plt
import timeit
import shutil
from scipy.stats import pearsonr
from scipy import stats
import scipy
import numpy
#disables x term  for use on atlas
from matplotlib.backends.backend_pdf import PdfPages
import json

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m','--mode',help='What is the mode you wish to run?')
    parser.add_argument('-pd','--pro_data',default='ProTherm.dat',help='The protherm .dat file')
    parser.add_argument('-d','--data',help='The data set to parse')
    parser.add_argument('-tp','--totalper',type=int,default=50,help='The total entries per mutation type desired in the benchmark set')
    parser.add_argument('-e','--entry',type=int,help='The entry in question')
    parser.add_argument('-bs','--benchmark_set',help='The file containing the keys for a previously picked benchmark set')
    parser.add_argument('-jn','--jobname',help='The name of the directory to place the job')
    parser.add_argument('-c','--cores',type=int,default=16,help='The number of cores to use in the ddg runs')
    parser.add_argument('-p','--pdbs',nargs="+",default=[],help='The list of pdb files. Names must contain the pdbid of the desired benchmark target')
    parser.add_argument('-rs','--runscript',default='runddg.sh',help='The submission script')
    parser.add_argument('-o','--output',help='The directory to output the prepped pdb files')
    parser.add_argument('-sk','--skipped',nargs="+",default=[],help='A way to manually skip certain entries based on pdbid to their quality or other concerns, currently only skips in the mut to pdbid mapping function')
    parser.add_argument('-re','--relaxex',default='/mnt/tank/bfrenz/home/Rosetta/main/source/bin/relax.default.linuxgccrelease',help='The path to the relax exectuable')
    parser.add_argument('-db','--database',default='/mnt/tank/bfrenz/home/Rosetta/main/database',help='The path to the Rosetta database')
    parser.add_argument('-sr','--skiprelax',type=int,default=0,help='Do you want to skip relax in the benchmark setup. ddG Results will be poor unless this is done manually later')
    parser.add_argument('-ms','--max_structure',type=int,default=50,help='The maximum number of mutations to choose that come from a single pdb file')
    parser.add_argument('-pt','--plottype',default='exact',help='What type of plot to make from the data')
    parser.add_argument('-t','--target_mutations',default='targetmutations.txt',help='The target binding ddg mutations')
    parser.add_argument('-sf','--scalefactor',type=float,default=1,help='What is the scale factor to apply to the scoring')
    parser.add_argument('-mut','--mutation',help='The mutation to search for')
    parser.add_argument('-ddg','--ddg',type=float,help='The ddg of the mutation to search for')
    parser.add_argument('-n','--num',help='The residue number to search for')
    parser.add_argument('-pmids','--pmids',nargs="+",help='The pmids to that require swapping')
    parser.add_argument('-js','--json',type=int,default=1,help='Are you working in json format? 1 for true')
    parser.add_argument('-a','--authors',nargs="+",default=[],help='The authors to search for')
    parser.add_argument('-oc','--outlier_cutoff',type=float,default='0.025',help='The change in pearsonsr required for a point to be considered an outlier')
    parser.add_argument('-mc','--matchcutoff',type=float,default=60,help='The cutoff required for two pdb sequences to be considered identical')
    parser.add_argument('-dups','--dups',help='The json file that matches pdbids to their duplicates')
    #binding ddg options
    parser.add_argument('-rw','--residue_window',type=int,default=4,help='The window of residues to extract')
    parser.add_argument('-ds','--dumpscript',default='dumpmuts.sh',help='The script that will dump the mutant pdb')
    parser.add_argument('-ss','--scorescript',default='score.sh',help='The script used to score the structures')
    parser.add_argument('-i','--index',type=int,help='The index in the json benchmarkset')
    parser.add_argument('-r','--relax',type=int,default=0,help='Set to 1 to relax the mutant pdb')
    #for dual job correlation
    parser.add_argument('-jn1','--job1',help='The directory of the first job')
    parser.add_argument('-jn2','--job2',help='The directory of the second job')
    args = parser.parse_args()
    return args

def parse_liz_set(args):
    mut = MUTATION_TYPES()
    stats = STATISTIC_TRACKING(mut)
    with open(args.data,'r') as liz_data:
        for line in liz_data:
            wt = line.split()[1]
            mut = line.split()[3]
            stats.add(wt,mut)
    print_stats(stats)

def parse_liz_set(args):
    mut = MUTATION_TYPES()
    stats = STATISTIC_TRACKING(mut)
    with open(args.data,'r') as liz_data:
        for line in liz_data:
            wt = line.split()[1]
            mut = line.split()[3]
            stats.add(wt,mut)
    print_stats(stats)

def parse_potapov(args):
    mut = MUTATION_TYPES()
    stats = STATISTIC_TRACKING(mut)
    num_muts = 0
    with open(args.data,'r') as potapov_data:
        for line in potapov_data:
            if line.startswith("#") or line.strip() == "":
                continue
            line = line.split(',')
            mut = line[2]
            muts = mut.split('_')
            for mut in muts:
                num_muts +=1
                mut = mut.split()
                wt = mut[1]
                mut = mut[3]
                stats.add(wt,mut)
    print 'total mutations,',num_muts
    print_stats(stats)

def parse_jia(args):
    mut = MUTATION_TYPES()
    stats = STATISTIC_TRACKING(mut)
    num_muts = 0
    with open(args.data,'r') as leidata:
        for line in leidata:
            if line.startswith('PDB_wild') or line.strip() == '':
                continue
            line = line.strip().split(',')
            wt = line[1]
            mut = line[3]
            stats.add(wt,mut)
            num_muts+=1
    print 'total mutations,',num_muts
    print_stats(stats)

def benchmark_stats(args):
    if args.json != 1:
        keys = parse_benchmark(args)
        protherm = parse_protherm(args)
        mut = MUTATION_TYPES()
        stats = STATISTIC_TRACKING(mut)
        for key in keys:
            wt = protherm[key].wild
            mut = protherm[key].mutation
            stats.add(wt,mut)
        print 'total mutations,',len(keys)
        print_stats(stats)
    else:
        benchmark = parse_benchmark_json(args)
        json_stats(args,benchmark)


#this is being depricated
def read_and_report_protherm_stats(args):
    
    mut = MUTATION_TYPES()
    stats = STATISTIC_TRACKING(mut)
    
    protherm_lines = ''
    with open(args.pro_data,'r') as pterm_file:
        protherm_lines= pterm_file.readlines()
    wt = ''
    mt = ''
    for line in protherm_lines:
        if line.startswith("MUTATION"):
            lline = line.split()
            if len(lline) < 2:
                continue
            if 'wild' in lline[1]:
                continue
            wt = lline[1]
            if len(lline) < 4:
                print lline
            mut = lline[3]
        if line.startswith('NO.'):
            wt = ''
            mut = ''
        if line.startswith('ddG'):
            if len(line.split()) < 2:
                continue
            ddg = line.split()[1]
            stats.add(wt,mut)
    print_stats(stats)

def parse_protherm(args):

    protherm_database = {}
    with open(args.pro_data,'r') as ptherm_file:
        lines = []
        modelid = 0
        for line in ptherm_file:
            if line.startswith('NO.'):
                if len(lines) < 1:
                    modelid = int(line.split()[1])
                    lines.append(line)
                    continue
                pdata = PROTHERM_ENTRY(modelid,lines)
                pdata.parse_lines()
                protherm_database[modelid] = pdata
                modelid = int(line.split()[1])
                lines = []
            lines.append(line)

    return protherm_database

def search_protherm(args):
    protherm = parse_protherm(args)
    #charged = ['K','R','D','E','H']
    for key in protherm.keys():
        entry = protherm[key]
        #if 'Campos' in " ".join(protherm[key].lines) and protherm[key].mutation in charged and entry.ddG != '':
        #    print " ".join(entry.lines)
        if entry.pdbwild == args.pdbs[0] and args.mutation == entry.mutation and args.num == entry.resnum and (entry.ddG == args.ddg or args.ddg == -99):
            print " ".join(entry.lines)
    
def parse_benchmark(args):
    keys = []
    with open(args.benchmark_set) as benchmark:
        for line in benchmark.readlines():
            if line.startswith("Entry Keys:"):
                line = line.split()
                it = 2
                while it < len(line):
                    keys.append(int(line[it]))
                    it+=1
    return keys

def parse_benchmark_json(args):
    pt = open(args.benchmark_set,'r')
    protherm = json.load(pt)
    pt.close()
    return protherm

def remove_unmatched_entries(args,benchmark):
    newbenchmark = deepcopy(benchmark)
    #for i in range(0,len(benchmark['data'])):
    newbenchmark['data'] = []
    for entry in benchmark['data']:
        pdbid = entry['PDBFileID']
        matched = match_pdb_to_list(args,pdbid)
        if matched != None:
            newbenchmark['data'].append(entry)
    print json.dumps(newbenchmark,sort_keys=True,indent=4,separators=(',',': '))
    
def setup_benchmark(args):
    # load data
    if not os.path.isdir(args.database) or not os.path.isfile(args.relaxex) and args.skiprelax == 0:
        print 'Either the database, the relax executable, or both are missing, please specify their locations with the flags -d and -re respectively. Alternatively skip relax with -sr'
    if args.json != 1:
        protherm = parse_protherm (args)
        keys = parse_benchmark(args)
        download_benchmark(args,protherm,keys)
        pdbids = get_pdbids(keys,protherm)
        mut_map = map_pdbids_to_mutations(args,pdbids,keys,protherm)
        #figure out target chains
        mut_map = determine_unassigned_chains(args,mut_map)
    else:
        benchmark = parse_benchmark_json(args)
        download_benchmark_json(benchmark)
        mut_map = json_to_mutmap(args,benchmark)
    
    
    pdbs = []
    if len(args.pdbs) > 0:
        pdbs = args.pdbs
    else:
        pdbs = glob.glob('pdbs/*.ent')
    args.pdbs = pdbs

    pdbnames = []
    if args.json == 1:
        pdbnames = clean_pdbs(args,benchmark)
    else:
        pdbnames = clean_pdbs(args,mut_map)
    print pdbnames
    
    relax_pdbs(args,pdbnames)

def clean_pdbs(args,benchmark):
    #create pdbs
    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    pdbnames = []
    if not args.json == 1:
        for pdbid in benchmark.keys():
            residues = pdbtools.get_unopened_residue_list(match_pdb_to_list(args,pdbid))
            pdbtools.convert_to_standard_aas(residues)
            residues = pdbtools.strip_non_protein(residues)
            for mutation in benchmark[pdbid]:
                dirname = args.output+'/'+pdbid+'_'+mutation.chain
                if not os.path.isdir(dirname):
                    os.mkdir(dirname)
                name = dirname+"/"+pdbid+"_"+mutation.chain+".pdb"
                if name not in pdbnames:
                    pdbnames.append((name))
                if not os.path.isfile(name):
                    chain_resis = pdbtools.get_chain_resis(residues,mutation.chain)
                    pdbtools.write_resis_to_pdb(chain_resis,name)
    else:
        finished_pdbs = []
        for entry in benchmark['data']:
            pdbid = entry['PDBFileID']
            #make sure we don't repeat pdbs we've already done
            if pdbid in finished_pdbs:
                continue
            residues = pdbtools.get_unopened_residue_list(match_pdb_to_list(args,pdbid))
            pdbtools.convert_to_standard_aas(residues)
            residues = pdbtools.strip_non_protein(residues)
            dirname = args.output+"/"+pdbid+"_"+entry['Mutations'][0]['Chain']
            name = dirname+'/'+pdbid+"_"+entry['Mutations'][0]['Chain']+".pdb"
            if not os.path.isdir(dirname):
                os.mkdir(dirname)
            if name not in pdbnames:
                pdbnames.append(name)
            if not os.path.isfile(name):
                chain_resis = pdbtools.get_chain_resis(residues,entry['Mutations'][0]['Chain'])
                pdbtools.write_resis_to_pdb(chain_resis,name)
    return pdbnames

def setup_pdb_directories(args):
    for pdb in args.pdbs:
        dirname = args.output+"/"+pdbtools.get_pdb_id(pdb)
        os.system('mkdir -p %s'%dirname)
        shutil.copy(pdb,dirname)

def relax_pdbs(args,pdbnames):

    if __name__ == '__main__':
        jobs = []
        for pdbname in pdbnames:
            jobs.append((args,pdbname))
        p = multiprocessing.Pool(args.cores,maxtasksperchild=1)
        p.map(que_relax, jobs)

def que_relax(inputs):
    args = inputs[0]
    pdbname = inputs[1]
    pdbid = pdbtools.get_pdb_id(pdbname)
    jobpath = args.output+"/"+pdbid
    shutil.copy('runrelax.sh',jobpath)
    shutil.copy('relax.xml',jobpath)
    self = os.path.realpath(__file__)
    command = ('nice -19 python %s -m "run relax" -p %s -o %s'%(self,pdbid,jobpath))
    os.system(command)

def relax_pdb(args):
    os.chdir(args.output)
    pdb = args.pdbs[0]+'.pdb'
    command = 'nice -19 sh runrelax.sh %s'%(pdb)
    os.system(command)
    pdbid = get_best_pdb()#pdbtools.get_pdb_id(pdbname)
    mvcommand = 'mv %s.pdb %s'%(pdbid,pdb)
    os.system('mv %s.pdb %s'%(pdbid,pdb))
    shutil.copy(pdb,'../')

def get_best_pdb():
    bestpdb = ''
    bestscore = 1e7
    with open('score.sc','r') as scorefile:
        for line in scorefile:
            if line.startswith("SCORE:"):
                line = line.split()
                if line[-1] == 'description':
                    continue
                score = float(line[1])
                if score < bestscore:
                    bestscore = score
                    bestpdb = line[-1]
    return bestpdb


def match_pdb_to_list(args,pdbid,chain=None):
    for pdb in args.pdbs:
        if chain == None:
            if pdbid.lower() in pdb.lower():
                return pdb
        else:
            to_match = pdbid.lower()+"_"+chain.lower()+".pdb"
            if to_match.lower() in pdb.lower():
                return pdb
    return None

def setup_pdbs(args,mut_map):
    for key in mut_map.keys():
        pdbs = glob.glob('pdbs/*')
        target = ''
        for pdb in pdbs:
            if key.lower in pdb:
                target = pdb
        resis = pdbtools.get_unopened_residue_list(target)
        multichain = is_multichain(resis)
        
def download_benchmark(args,protherm,keys):
    pdbids = []
    for key in keys:
        if protherm[key].pdbwild == ' ':
            print key
        pdbids.append(protherm[key].pdbwild)
    for pdbid in pdbids:
        download_pdb(pdbid)
    return protherm

def download_benchmark_json(jsonbenchmark):
    for i in range(0,len(jsonbenchmark['data'])-1):
        pdbid = jsonbenchmark['data'][i]['PDBFileID']
        download_pdb(pdbid)
    return

def get_pdbids(keys,protherm):
    pdbids = []
    for key in keys:
        if protherm[key].pdbwild in pdbids:
            continue
        pdbids.append(protherm[key].pdbwild)
    return pdbids

def run_benchmark(args):
    mut_map = {}
    if args.json != 1:
        keys = parse_benchmark(args)
        protherm = parse_protherm(args)
        pdbids = get_pdbids(keys,protherm)
        mut_map = map_pdbids_to_mutations(args,pdbids,keys,protherm)
        mut_map = determine_unassigned_chains(args,mut_map)
    else:
        benchmark = parse_benchmark_json(args)
        mut_map = json_to_mutmap(args,benchmark)
        verify_pdbs(args,benchmark)
    mut_map = convert_map_to_rosetta_num(args,mut_map)
    original_mut_map = mut_map
    ddg_data = load_ddg_data(args)
    mut_map = remove_completed_jobs(mut_map,ddg_data)
    
    #count the total mutations
    total_jobs = 0
    mutation_jobs = []
    for key in mut_map.keys():
        unique_chains = []
        for mut in mut_map[key]:
            if mut.chain not in unique_chains:
                total_jobs+=1
                mutation_jobs.append((key,mut.chain))
                unique_chains.append(mut.chain)
    print len(mutation_jobs)


    if __name__ == '__main__':
        jobcount = 0
        p = multiprocessing.Pool(args.cores)
        jobs = []
        while jobcount < total_jobs:
            pdbid = mutation_jobs[jobcount][0]
            chain = mutation_jobs[jobcount][1]
            pdb = match_pdb_to_list(args,pdbid,chain)
            jobs.append((args,pdb,mut_map))
            jobcount+=1
        p.map(submit_ddg, jobs)

def verify_pdbs(args,benchmark):
    all_present = True
    for entry in benchmark['data']:
        pdbname = entry['PDBFileID']+"_"+entry['Mutations'][0]['Chain']
        matched = match_pdb_to_list(args,pdbname)
        if matched == None:
            all_present = False
            print pdbname,'is missing exiting protocol'
    if not all_present:
        exit()
    return

def submit_ddg(inputs):
    args = inputs[0]
    pdb = inputs[1]
    mut_map = inputs[2]
    pdbid = pdbtools.get_pdb_id(pdb)
    cdir = os.getcwd()
    if not os.path.isdir(args.jobname):
        os.mkdir(args.jobname)
    if not os.path.isdir('%s/%s'%(args.jobname,pdbid)):
        os.mkdir('%s/%s'%(args.jobname,pdbid))
    path = "%s/%s/"%(args.jobname,pdbid)
    chain = pdbid.split("_")[1]
    mutmapkey = pdbid.split("_")[0]
    mutations = mut_map[mutmapkey]
    #skip jobs that are already finished
    if os.path.isfile('%smutfile.ddg'%path):
        print path,'is finished'
        return
    create_mutfile(mutations,chain,'%smutfile.mut'%path)
    self = os.path.realpath(__file__)
    command = 'nice -19 python %s -m pjob -p %s -jn %s -rs %s '%(self,pdb,args.jobname,args.runscript)
    print command
    os.system(command)

def create_mutfile(mutations,chain,name):
    with open(name,'w') as mutfile:
        total_for_chain = 0
        for mutation in mutations:
            if mutation.chain == chain:
                total_for_chain+=1
        mutfile.write('total %s'%total_for_chain+"\n")
        #mutfile.write(str(len(mutations))+"\n")
        for mutation in mutations:
            if mutation.chain != chain:
                continue
            mutfile.write(str(1)+"\n")
            mutfile.write("%s %s %s"%(mutation.wt,mutation.resnum,mutation.mut)+"\n")
    
def run_ddg(args):
    target_pdb = args.pdbs[0]
    pdbid = pdbtools.get_pdb_id(target_pdb)
    
    os.system('cp %s %s %s/%s/'%(args.runscript,target_pdb,args.jobname,pdbid))
    os.chdir("%s/%s"%(args.jobname,pdbid))
    target_pdb = target_pdb.split('/')[-1]
    
    if target_pdb == '':
       print 'no proper pdbfile found'
       exit()
    command = 'nice -19 sh %s %s'%(args.runscript,target_pdb)
    os.system(command)

def map_pdbids_to_mutations(args,pdbids,keys,protherm):
    pdbid_map = defaultdict(list)
    for key in keys:
        mutation = Mutation(protherm[key].resnum,protherm[key].mutchain,protherm[key].wild,protherm[key].mutation)
        mutation.set_protherm_key(key)
        if protherm[key].pdbwild in args.skipped:
            continue
        pdbid_map[protherm[key].pdbwild].append(mutation)
    return pdbid_map

def json_to_mutmap(args,jsonbenchmark):
    pdbid_map = defaultdict(list)
    for i in range(0,len(jsonbenchmark['data'])-1):
        datapoint = jsonbenchmark['data'][i]['Mutations'][0]
        mutation = Mutation(datapoint['ResidueID'],datapoint['Chain'],datapoint['WildTypeAA'],datapoint['MutantAA'])
        mutation.set_protherm_key(i)
        if jsonbenchmark['data'][i]['PDBFileID'] in args.skipped:
            continue
        pdbid_map[ jsonbenchmark['data'][i]['PDBFileID'] ].append(mutation)
    return pdbid_map

def load_ddg_data(args,jobname=None):
    
    if jobname == None:
        jobname = args.jobname
    if jobname == None:
        print 'error attempting to load ddg data and no job name is specificied. Please use the flag -jn to set the job name to the parent directory containing each of your jobs'
        exit()
    #read the outfile and collect the average scores for each mutation
    ddg_data = defaultdict(dict)
    for job in glob.glob(jobname+"*/"):
        pdb_ddgs = defaultdict(list)
        mutfile = job+'mutfile.ddg'
        if not os.path.isfile(mutfile):
            continue
        with open(mutfile) as mfile:
            for line in mfile:
                data = line.split()
                mutation = data[2].replace(':','')
                if mutation != 'WT':
                    mutation = mutation.split('_')[1]
                    num = int(re.split('(\d+)',mutation)[1])
                    residue = re.split('(\d+)',mutation)[2]
                    residue = amino_acids.longer_names[residue]
                    score = float(data[3])
                    pdb_ddgs[(num,residue)].append(score)
                if mutation == 'WT':
                    score = float(data[3])
                    pdb_ddgs[mutation].append(score)
        pdbid = re.split('/|_',job)[-3]
        chain = re.split('/|_',job)[-2]
        ddg_data[(pdbid,chain)] = pdb_ddgs



    for key in ddg_data.keys():
        for mutation in ddg_data[key].keys():
            average_score = sum(ddg_data[key][mutation])/len(ddg_data[key][mutation])
            if len(ddg_data[key][mutation]) < 3:
                print 'WARNING: fewer than 3 data points for ',key,mutation,'removing it from the ddg list'
                del ddg_data[key][mutation]
            else:
                ddg_data[key][mutation] = average_score

    #loop over the loaded data and convert average scores to ddgs
    for pdbid in ddg_data.keys():
        if 'WT' not in ddg_data[pdbid].keys():
            print 'Warning there is no wild type score for',pdbid
            exit()
        for mutation in ddg_data[pdbid].keys():
            if mutation == 'WT':
                continue
            rosddg = (ddg_data[pdbid][mutation]-ddg_data[pdbid]['WT'])/args.scalefactor
            ddg_data[pdbid][mutation] = rosddg

    #print ddgs
    #for pdbid in ddg_data.keys():
    #    for mutation in ddg_data[pdbid].keys():
    #        print pdbid,mutation,ddg_data[pdbid][mutation]
    return ddg_data

def load_p3_data(args):
    if args.jobname == None:
        print 'error attempting to load ddg data and no job name is specificied. Please use the flag -jn to set the job name to the parent directory containing each of your jobs'
        exit()
    #read the outfile and collect the average scores for each mutation
    ddg_data = defaultdict(dict)
    for job in glob.glob(args.jobname+"*/"):
        pdb_ddgs = defaultdict(list)
        predictions = job+'ddg_predictions.out'
        if not os.path.isfile(predictions):
            continue
        with open(predictions) as pfile:
            for line in pfile:
                if 'description' in line or line.strip() == '':
                    continue
                data = line.split()
                target = data[1]
                target = re.split('(\D+)',target)
                wild = target[1]
                num = int(target[2])
                mut = target[3]
                score = float(data[2])
                pdb_ddgs[(num,mut)] = score
        pdbid = re.split('/|_',job)[-3]
        chain = re.split('/|_',job)[-2]
        ddg_data[(pdbid,chain)] = pdb_ddgs
    return ddg_data

#this will assign a chain id to every mutation. If the model is a homodimer the first chain is selected if the model contains multiple chains and no chain is
#specificed in protherm it will assign the chain that matches the reported wild type sequence. If something other than a single chain is identified an erorr will be reported.
def determine_unassigned_chains(args,mut_map):
    new_mut_map = defaultdict(list)
    for pdbid in mut_map.keys():
        residues = pdbtools.get_unopened_residue_list(match_pdb_to_list(args,pdbid))
        pdbtools.convert_to_standard_aas(residues)
        residues = pdbtools.strip_non_protein(residues)
        chains = pdbtools.get_chains(residues)
        for mutation in mut_map[pdbid]:
            if mutation.chain in chains:
                new_mut_map[pdbid].append(mutation)
                continue
            if mutation.chain == '':
                if len(chains) > 1:
                    print "ERROR AMBIGUOUS CHAINS IN",pdbid
                else:
                    mutation.chain = chains[0]
            else:
                #this block of text is to find the non redundant sequences and then finds the chain that matches the native mutation
                sequences = pdbtools.get_sequences(residues)
                non_redundant = []
                for sequence in sequences:
                    if sequence not in non_redundant:
                        non_redundant.append(sequence)
                target_chains = []
                desired_chains = []
                for i in range(0,len(non_redundant)):
                    target_chains.append(chains[i])
                for res in residues:
                    if res.chain in target_chains and res.num == int(mutation.resnum) and amino_acids.longer_names[res.name] == mutation.wt:
                        desired_chains.append(res.chain)
                if len(desired_chains) > 1:
                    matched_chains = read_pdb_chains(match_pdb_to_list(args,pdbid))
                    unique_chains = []
                    identical_chains = []
                    for chain in desired_chains:
                        if chain not in identical_chains:
                            unique_chains.append(chain)
                            for mchains in matched_chains:
                                if chain in mchains:
                                    for mchain in mchains:
                                        if mchain not in identical_chains:
                                            identical_chains.append(mchain)
                    desired_chains = unique_chains
                if len(desired_chains) != 1:
                    print pdbid,mutation.resnum,mutation.wt,'Error: ',len(desired_chains),'found. It should have 1. Skipping'
                    continue
                mutation.chain = desired_chains[0]
                new_mut_map[pdbid].append(mutation)
    return new_mut_map

#converts the map so that Mutation.resnum is in the rosettanumbering. Stores the original value as mutation.fnum
def convert_map_to_rosetta_num(args,pdbid_map):
    new_map = defaultdict(list)
    for key in pdbid_map:
        for mutation in pdbid_map[key]:
            matched_pdb = match_pdb_to_list(args,key,mutation.chain)
            resis = pdbtools.get_unopened_residue_list(matched_pdb)
            counter = 1
            for res in resis:
                mutationid = re.split('(\D+)',mutation.resnum)
                mutationnum = mutationid[0]
                mutationicode = ' '
                if len(mutationid) > 1:
                    mutationicode = mutationid[1]
                if res.chain == mutation.chain and res.num == int(mutationnum) and mutationicode == res.icode:
                    newMutation = Mutation(counter,mutation.chain,mutation.wt,mutation.mut)
                    newMutation.set_pnum(mutation.resnum)
                    newMutation.set_protherm_key(mutation.get_protherm_key())
                    new_map[key].append(newMutation)
                counter+=1
    for key in pdbid_map.keys():
        if key not in new_map.keys():
            print key,pdbid_map[key][0].chain
    assert len(pdbid_map.keys()) == len(new_map.keys())
    for key in pdbid_map.keys():
        if len(pdbid_map[key]) != len(new_map[key]):
            print 'Warning: some mutations could not be found in for pdb',key,'The are being skipped'
            print len(pdbid_map[key]),len(new_map[key])
    return new_map

def add_rosettanum_to_json(args,benchmark):
    for entry in benchmark['data']:
        added_rosnum = False
        pdbid = entry['PDBFileID']
        chain = entry['Mutations'][0]['Chain']
        matched_pdb = match_pdb_to_list(args,pdbid,chain)
        resis = pdbtools.get_unopened_residue_list(matched_pdb)
        counter = 1
        mutationid = entry['Mutations'][0]['ResidueID']
        mutationnum,mutationicode = split_mutationid(mutationid)
        for res in resis:
            if res.chain == chain and res.num == int(mutationnum) and mutationicode == res.icode:
                entry['Mutations'][0]['RosettaNum'] = counter
                added_rosnum = True
            counter+=1
        if added_rosnum == False:
            print 'Could not find mutation for',pdbid,chain,mutationnum,mutationicode
            exit()
    return benchmark


def add_rosettanum_to_protherm(args):
    protherm = parse_protherm_json(args)
    for entry in protherm['data']:
        added_rosnum = False
        pdbid = entry['PDBFileID']
        pdbid = pdbid.lower()
        chain = entry['Mutations'][0]['Chain']
        for pdb in args.pdbs:
            if pdbid in pdb:
                matched_pdb = pdb
        resis = pdbtools.get_unopened_residue_list(matched_pdb)
        counter = 1
        mutationid = entry['Mutations'][0]['ResidueID']
        mutationnum,mutationicode = split_mutationid(mutationid)
        for res in resis:
            if res.chain == chain and res.num == int(mutationnum) and mutationicode == res.icode:
                entry['Mutations'][0]['RosettaNum'] = counter
                added_rosnum = True
            if res.chain  == chain:
                counter+=1
        if added_rosnum == False:
            print 'Could not find mutation for',pdbid,chain,mutationnum,mutationicode
            exit()
    print json.dumps(protherm,sort_keys=True,indent=4,separators=(',',': '))

def remove_completed_jobs(mutation_map,ddg_data):
    new_mut_map = defaultdict(list)
    for pdbid in mutation_map.keys():
        for mutation in mutation_map[pdbid]:
            completed = False
            if (pdbid,mutation.chain) in ddg_data.keys():
                if (mutation.resnum,mutation.mut) in ddg_data[(pdbid,mutation.chain)].keys():
                    completed = True
            if completed:
                print pdbid,mutation.chain,mutation.resnum,mutation.mut,'is finished'
            else:
                new_mut_map[pdbid].append(mutation)
    return new_mut_map

def split_mutationid(mutationid):
    mutationid = re.split('(\D+)',mutationid)
    mutationnum = mutationid[0]
    mutationicode = ' '
    if len(mutationid) > 1:
        mutationicode = mutationid[1]
    return mutationnum,mutationicode

def is_multichain(residues):
    chains = []
    for residue in residues:
        if residue.chain not in chains:
            chains.append(residue.chain)
        if len(chains) > 1:
            return True
    return False

def strip_monomer(residues):
    chains = []
    monomer_res = []
    for residue in residues:
        if residue.name not in amino_acids.longer_names:
            continue
        if residue.chain not in chains and len(chains) == 0:
            chains.append(residue.chain)
        if residue.chain in chains:
            monomer_res.append(residue)
    return monomer_res


#reads the chains that are associated with the same compound to figure out what is grouped
def read_pdb_chains(pdbname):
    matched_chains = []
    with open(pdbname,'r') as pdbfile:
        for line in pdbfile:
            if line.startswith('COMPND') and 'CHAIN:' in line:
                line = re.split('CHAIN:|;',line)
                chains = line[1]
                chains = chains.split()
                grouped_chains = []
                for chain in chains:
                    chain = chain.strip(',')
                    grouped_chains.append(chain)
                matched_chains.append(grouped_chains)
    return matched_chains


def print_stats(stats):
    print 'positive to negative,',stats.pos2neg
    print 'negative to positive,',stats.neg2pos
    print 'positive to polar non-charged,', stats.pos2ncp
    print 'negative to poloar non-charged,', stats.neg2ncp
    print 'polar non-charged to positive,', stats.ncp2pos
    print 'polar non-charged to negative,', stats.ncp2neg
    print 'negative to hydrophobic,', stats.neg2h
    print 'hydrophobic to negative,', stats.h2neg
    print 'positive to hydrophobic,', stats.pos2h
    print 'hydrophobic to positive,', stats.h2pos
    print 'non-charged polar to hydrophobic,', stats.ncp2h
    print 'hydrophobic to non-charged polar,', stats.h2ncp
    print 'non-charged polar to non-charged polar,', stats.ncp2ncp
    print 'hydrophobic to hydrophobic,', stats.h2h
    print 'like to like charge,', stats.crg2crg
    print 'involves proline,', stats.pro
    print 'involves cysteine,', stats.cys
    print 'small to large,',stats.sm2l
    print 'large to small,',stats.l2sm
    print 'same size to same size,', stats.siz2siz

def print_entry_lines(args):
    protherm = parse_protherm(args)
    print " ".join(protherm[args.entry].lines)
    print protherm[args.entry].wild,protherm[args.entry].mutation

def count_ddg_data(args):
    total_ddg = 0
    total_pdb = 0
    with open(args.pro_data,'r') as ptherm:
        haspdb = False
        for line in ptherm:
            if line.startswith('NO.'):
                haspdb = False
            if line.startswith("PDB_wild"):
                lline = line.split()
                if len(lline) > 1:
                    haspdb = True
            if line.startswith("ddG"):
                lline = line.split()
                if len(lline) > 1:
                    total_ddg+=1
                    if haspdb:
                        total_pdb+=1
    print 'there are',total_pdb,'entries with pdbs and ddg data. there are',total_ddg,'total entries with ddg data'

#returns the info excluding the first word of the line
def get_line_data(line):
    line = line.split()
    if len(line) > 1:
        return ' '.join(line[1:])
    else:
        return ' '

def download_pdb(pdbid):
    pdbl = PDBList()
    pdbfile = "pdb"+pdbid+".ent"
    if not os.path.isfile('pdbs/%s'%pdbfile):
        try:
            pdbl.retrieve_pdb_file(pdbid,pdir='pdbs')
        except:
            print 'could not download', pdbid

def filter_nonddG_entries(protherm):
    filtered_protherm = {}
    for key in protherm.keys():
        if protherm[key].ddG != '':
            filtered_protherm[key] = protherm[key]
    return filtered_protherm

#picks a benchmark set out of the protherm database
def pick_from_protherm(args):
    protherm = parse_protherm(args)
    protherm = filter_nonddG_entries(protherm)

    stats = STATISTIC_TRACKING(MUTATION_TYPES())
    done = False
    benchmarkset = []
    n = 100000
    wild_pdbs = []
    acceptance_level = 1 #1 will accept only unique categories. 2 will accept prolines and cysteines, 3 will accept any unmet broader categories
    while not done:
        newstats = deepcopy(stats)
        entry = random.choice(protherm.keys())
        n-=1
        if n < 0:
            done = True
        if n == 0:
            done = True
        if protherm[entry].pdbwild in args.skipped:
            continue
        
        #check if entry is already accepted and valid
        if entry in benchmarkset:
            continue
        if protherm[entry].mutation not in amino_acids.amino_acids or protherm[entry].wild not in amino_acids.amino_acids or protherm[entry].pdbwild == ' ':
            continue

        if wild_pdbs.count(protherm[entry].pdbwild) > args.max_structure:
            continue

        
        newstats.add(protherm[entry])
        if newstats.added_new(stats,args.totalper,acceptance_level):
            if not is_duplicate(benchmarkset,entry,protherm):
                stats = deepcopy(newstats)
                wild_pdbs.append(protherm[entry].pdbwild)
                benchmarkset.append(entry)
        if stats.has_enough(args.totalper,acceptance_level):
            if acceptance_level < 3:
                acceptance_level+=1
            else:
                done = True

    #after searching for random options loop over the entire protherm to fill the missing categories
    acceptance_level = 1
    while(acceptance_level <= 3):
        for entry in protherm.keys():
            newstats = deepcopy(stats)
            if protherm[entry].pdbwild in args.skipped:
                continue
            
            #check if entry is already accepted and valid
            if entry in benchmarkset:
                continue
            if protherm[entry].mutation not in amino_acids.amino_acids or protherm[entry].wild not in amino_acids.amino_acids or protherm[entry].pdbwild == ' ':
                continue

            #if wild_pdbs.count(protherm[entry].pdbwild) > args.max_structure:
            #    continue

            newstats.add(protherm[entry])
            if newstats.added_new(stats,args.totalper,acceptance_level):
                if not is_duplicate(benchmarkset,entry,protherm):
                    stats = deepcopy(newstats)
                    wild_pdbs.append(protherm[entry].pdbwild)
                    benchmarkset.append(entry)
            if stats.has_enough(args.totalper,acceptance_level):
                acceptance_level +=1
                break
        acceptance_level+=1
    print_stats(stats)
    print 'Entry Keys:',' '.join(str(x) for x in benchmarkset)

#picks a benchmark set out of the protherm database
def pick_from_curated_protherm(args):
    protherm = parse_protherm_json(args)
    duplicates = {}
    if args.dups != None:
        dp = open(args.dups,'r')
        duplicates = json.load(dp)
        dp.close()
    if args.output != None:
        download_benchmark_json(protherm)
    benchmark = {}
    benchmark['data'] = list({})
    benchmark = clean_json(benchmark)

    stats = STATISTIC_TRACKING(MUTATION_TYPES())
    done = False
    benchmarkset = []
    n = 100000
    wild_pdbs = []
    acceptance_level = 1 #1 will accept only unique categories. 2 will accept prolines and cysteines, 3 will accept any unmet broader categories
    while not done:
        newstats = deepcopy(stats)
        entry = random.randint(0,len(protherm['data'])-1)
        n-=1
        if n < 0:
            done = True
        if n == 0:
            done = True
        if protherm['data'][entry]['PDBFileID'] in args.skipped:
            continue
        #check if entry is already accepted and valid
        pdbid = protherm['data'][entry]['PDBFileID'].lower()
        chain = protherm['data'][entry]['Mutations'][0]['Chain']
        resid = protherm['data'][entry]['Mutations'][0]['ResidueID']
        mutantAA = protherm['data'][entry]['Mutations'][0]['MutantAA']
        try:
            float(resid)
            resid=resid+' '
        except:
            sallgood = True
        rosnum = protherm['data'][entry]['Mutations'][0]['RosettaNum']
        matches = find_matching_mutations(duplicates,pdbid,chain,resid,mutantAA)
        matches.append((pdbid,chain,rosnum,mutantAA))
        skip = False
        for match in matches:
            if match in benchmarkset:
                skip  = True
                break
        if skip:
            continue

        if wild_pdbs.count(protherm['data'][entry]['PDBFileID']) > args.max_structure:
            continue
        
        newstats.add_json(protherm['data'][entry]['Mutations'][0])
        if newstats.added_new(stats,args.totalper,acceptance_level):
            stats = deepcopy(newstats)
            wild_pdbs.append(protherm['data'][entry]['PDBFileID'])
            best_experiment = pick_best_experiment(protherm['data'][entry]['ExperimentalDDGs'])
            benchmark['data'].append(protherm['data'][entry])
            benchmark['data'][-1]['ExperimentalDDGs'] = best_experiment
            for match in matches:
                benchmarkset.append(match)
        if stats.has_enough(args.totalper,acceptance_level):
            if acceptance_level < 3:
                acceptance_level+=1
            else:
                done = True

    #after searching for random options loop over the entire protherm to fill the missing categories
    acceptance_level = 1
    print len(benchmark['data']),'entries after first pass'
    while(acceptance_level <= 3):
        for i in range(0,len(protherm['data'])):
            newstats = deepcopy(stats)
            if protherm['data'][i]['PDBFileID'] in args.skipped:
                continue
            
            #check if entry is already accepted and valid
            pdbid = protherm['data'][i]['PDBFileID'].lower()
            chain = protherm['data'][i]['Mutations'][0]['Chain']
            resid = protherm['data'][i]['Mutations'][0]['ResidueID']
            mutantAA = protherm['data'][i]['Mutations'][0]['MutantAA']
            try:
                float(resid)
                resid=resid+' '
            except:
                sallgood = True
            rosnum = protherm['data'][i]['Mutations'][0]['RosettaNum']
            matches = find_matching_mutations(duplicates,pdbid,chain,resid,mutantAA)
            matches.append((pdbid,chain,rosnum,mutantAA))
            skip = False
            for match in matches:
                if match in benchmarkset:
                    skip = True
                    break
            if skip:
                continue
            
            newstats.add_json(protherm['data'][i]['Mutations'][0])
            if newstats.added_new(stats,args.totalper,acceptance_level):
                stats = deepcopy(newstats)
                best_experiment = pick_best_experiment(protherm['data'][i]['ExperimentalDDGs'])
                benchmark['data'].append(protherm['data'][i])
                benchmark['data'][-1]['ExperimentalDDGs'] = best_experiment
                for match in matches:
                    benchmarkset.append(match)
            if stats.has_enough(args.totalper,acceptance_level):
                acceptance_level +=1
        acceptance_level+=1
    benchmark['stats'] = stats.to_dict()
    #print_stats(stats)
    #download_benchmark_json(benchmark)
    pdbs = []
    if len(args.pdbs) > 0:
        pdbs = args.pdbs
    else:
        pdbs = glob.glob('pdbs/*.ent')
    args.pdbs = pdbs
    args.pdbs = clean_pdbs(args,benchmark)
    print len(benchmark['data']),'entries total'

    with open(args.benchmark_set,'w') as newbenchmark:
        newbenchmark.write(json.dumps(benchmark,sort_keys=True,indent=4,separators=(',',': ')))

def load_pdbs(args):
    pdbmap = defaultdict(dict)
    for pdb in args.pdbs:
        resis = pdbtools.get_unopened_residue_list(pdb)
        pdbid = re.split('pdb|.ent',pdb)[-2]
        sequences = pdbtools.get_sequences(resis)
        pdbmap[pdbid]['sequences'] = []
        chains = []
        pdbmap[pdbid]['residues'] = defaultdict(list)
        for res in resis:
            if res.chain not in chains:
                chains.append(res.chain)
            pdbmap[pdbid]['residues'][res.chain].append(res)
        for i in range(0,len(chains)):
            xs = 0.0
            sequence = sequences[i]
            aaseq = []
            for aa in sequence:
                if aa not in amino_acids.amino_acids:
                    continue
                aaseq.append(aa)
            if len(aaseq) == 0:
                continue
            aaseq = ''.join(aaseq)
            pdbmap[pdbid]['sequences'].append((chains[i],aaseq))
    return pdbmap

def find_duplicates(args):
    pdbs = load_pdbs(args)
    duplicates = defaultdict(list)# defaultdict(defaultdict(defaultdict(list)))
    for i in range(0,len(pdbs.keys())):
        key = sorted(pdbs.keys())[i]
        if key not in duplicates.keys():
            duplicates[key] = {}
        for match1 in pdbs[key]['sequences']:
            chain = match1[0]
            seq1 = match1[1]
            try:
                len(duplicates[key][chain])
            except:
                duplicates[key][chain] = []
            for ii in range(i+1,len(pdbs.keys())):
                key2 = sorted(pdbs.keys())[ii]
                if key2 not in duplicates.keys():
                    duplicates[key2] = {}
                for match2 in pdbs[key2]['sequences']:
                    chain2 = match2[0]
                    try:
                        len(duplicates[key2][chain2])
                    except:
                        duplicates[key2][chain2] = []
                    seq2 = match2[1]
                    duplicate,string1,string2 = is_duplicate_sequence(args,seq1,seq2)
                    if duplicate:
                        print duplicate
                        num_alignment1 = match_res_to_alignment(pdbs[key]['residues'][chain],string1)
                        num_alignment2 = match_res_to_alignment(pdbs[key2]['residues'][chain2],string2)
                        duplicates[key][chain].append({'chain':chain2,'pdbid':key2,'alignment':num_alignment1})
                        duplicates[key2][chain2].append({'chain':chain,'pdbid':key,'alignment':num_alignment2})
    with open(args.output,'w') as outfile:
        outfile.write(json.dumps(duplicates,sort_keys=True,indent=4,separators=(',',': ')))

def find_matching_mutations(duplicates,pdb1,chain1,resid,mutantAA):
    matches = []
    if pdb1 not in duplicates.keys():
        return matches
    for duplicate in duplicates[pdb1][chain1]:
        matched_pdb = duplicate['pdbid']
        matched_chain = duplicate['chain']
        alignment_res = duplicate['alignment'][resid]
        matches.append((matched_pdb,matched_chain,alignment_res,mutantAA))
    return matches

def remove_interfaces(args):
    protherm = parse_protherm_json(args)
    newprotherm = {'data':[]}
    clean_json(newprotherm)
    pdbs = load_pdbs(args)
    for entry in protherm['data']:
        for mutation in entry['Mutations']:
            chain = mutation['Chain']
            try :
                resnum = int(mutation['ResidueID'])
                icode = ' '
            except:
                resnum = int(mutation['ResidueID'][:-1])
                icode = mutation['ResidueID'][-1]
            pdbid = entry['PDBFileID'].lower()
            matched_res = None
            for chainkey in pdbs[pdbid]['residues'].keys():
                if chainkey != chain:
                    continue
                for res in pdbs[pdbid]['residues'][chainkey]:
                    if res.chain == chain and resnum == res.num and icode == res.icode:
                        matched_res = res
                        break
            if not is_interface(matched_res,pdbs[pdbid]['residues']):
                newprotherm['data'].append(entry)
    print json.dumps(newprotherm,sort_keys=True,indent=4,separators=(',',': '))

def match_res_to_alignment(residues,alignstring):
    resid_to_alignment = {}
    rescounter = 0
    for i in range(0,len(alignstring)):
        if not alignstring[i] == '-':
            resid = str(residues[rescounter].num)+residues[rescounter].icode
            resid_to_alignment[resid] = i
            rescounter+=1
    return resid_to_alignment

def is_interface(res,chains):
    for chain in chains:
        if chain == res.chain:
            continue
        for matchedres in chains[chain]:
            for atom1 in res.atoms:
                for atom2 in matchedres.atoms:
                    if pdbtools.atom_dist(atom1,atom2) < 5:
                        return True
    return False

def in_benchmark(dups,benchmark,pdbid,chain,resid):
    for entry in benchmark:
        bm_pdbid = entry['PDBFileID']
        
    #for dup in dups[pdbid][chain]:
    #    if 
    


def is_duplicate_sequence(args,seq1,seq2):
    alignments = pairwise2.align.localms(seq1,seq2,2,1,-.5,-.1)
    alignseq1 = alignments[0][0]
    alignseq2 = alignments[0][1]
    matching = 0.0
    for i in range(0,len(alignseq1)):
        if alignseq1[i] == alignseq2[i] and alignseq1[i] != '-':
            matching+=1
    percentMatching = matching/len(alignseq1)*100
    if percentMatching > args.matchcutoff:
        return True,alignseq1,alignseq2
    else:
        return False
    

def parse_protherm_json(args):
    pt = open(args.pro_data,'r')
    protherm = json.load(pt)
    pt.close()
    return protherm

def json_stats(args,protherm):
    duplicates = {}
    if args.dups != None:
        dp = open(args.dups,'r')
        duplicates = json.load(dp)
        dp.close()
    uniques = []
    mut = MUTATION_TYPES()
    stats = STATISTIC_TRACKING(mut)
    repeats = []
    for i in range(0,len(protherm['data'])-1):
        mutation = protherm['data'][i]['Mutations'][0]
        pdbid = protherm['data'][i]['PDBFileID'].lower()
        chain = mutation['Chain']
        resid = mutation['ResidueID']
        mutantAA = mutation['MutantAA']
        try:
            float(resid)
            resid=resid+' '
        except:
            sallgood = True
        rosnum = mutation['RosettaNum']
        matches = find_matching_mutations(duplicates,pdbid,chain,resid,mutantAA)
        matches.append((pdbid,chain,rosnum,mutantAA))
        skip = False
        for match in matches:
            if match in repeats:
                skip = True
                break
        if skip:
            print 'skipping'
            continue
        if( (protherm['data'][i]['PDBFileID'],mutation['MutantAA'],mutation['WildTypeAA'],mutation['Chain'],mutation['ResidueID']) not in uniques):
            uniques.append((protherm['data'][i]['PDBFileID'],mutation['MutantAA'],mutation['WildTypeAA'],mutation['Chain'],mutation['ResidueID']))
            stats.add_json(mutation)
        for match in matches:
            repeats.append(match)

    print len(uniques)
    print_stats(stats)

def clean_json(benchmark):
    for x in benchmark['data']:
        benchmark['data'][x]['ExperimentalDDGs'] = {}
    return benchmark

def pick_best_experiment(experiments):
    best_value = 1e5
    best_experiment = {}
    if not isinstance(experiments,list):
        return experiments
    for experiment in experiments:
        ph = float(experiment['pH'])
        if abs(7-ph) < best_value:
            best_value = abs(7-ph)
            best_experiment = experiment
    return best_experiment 

def is_duplicate(benchmark,key,protherm):
    pdb = protherm[key].pdbwild
    residue = protherm[key].resnum
    mutation = protherm[key].mutation
    for element in benchmark:
        bpro = protherm[element]
        if bpro.pdbwild == pdb and bpro.resnum == residue and bpro.mutation == mutation:
            return True
    return False

def protherm_stats(args):
    protherm = parse_protherm(args)
    protherm = filter_nonddG_entries(protherm)
    #print len(protherm.keys()),'entries after filter'

    stats = STATISTIC_TRACKING( MUTATION_TYPES() ) 
    uncounted_entries = 0
    counted_entries = 0
    multi_mutation = 0
    nonduplicates = []
    total_muts = 0
    for entry in protherm:
        if is_duplicate(nonduplicates,entry,protherm):
            continue
        if protherm[entry].pdbwild == ' ':
            continue
        total_muts +=1
        newstats = deepcopy(stats)
        newstats.add(protherm[entry])
        nonduplicates.append(entry)
        if protherm[entry].mutation == 'multi':
            multi_mutation+=1
        if not newstats.added_new(stats,100000):
            #if protherm[entry].mutation != 'multi' and protherm[entry].wild != 'wild':
            uncounted_entries+=1
        else:
            counted_entries+=1
        stats = newstats
    print 'total mutations,',total_muts
    print_stats(stats)
    #print uncounted_entries,'are not included in the statistics.',counted_entries,'are',multi_mutation,'contain multiple mutations'

def validate_mutations(stats):
    for aa1 in amino_acids.amino_acids:
        for aa2 in amino_acids.amino_acids:
            newstats = deepcopy(stats)
            newstats.add(aa1,aa2)
            if not newstats.added_new(stats,1e10):
                #print_stats(stats)
                #print ' '
                #print_stats(newstats)
                print aa1,aa2,'was not added'

def plot_ddg_data(args,p3=False):
    print 'loading data...'
    ddg_data = {}
    if p3 == False:
        ddg_data = load_ddg_data(args)
    else:
        ddg_data = load_p3_data(args)

    protherm_ddgs = []
    rosetta_ddgs = []
    results = []
    mut_map = {}
    if not args.json == 1:
        keys = parse_benchmark(args)
        protherm = parse_protherm(args)
        pdbids = get_pdbids(keys,protherm)
        mut_map = map_pdbids_to_mutations(args,pdbids,keys,protherm)
        mut_map = determine_unassigned_chains(args,mut_map)
        mut_map = convert_map_to_rosetta_num(args,mut_map)

        for key in keys:
            if  protherm[key].pdbwild not in mut_map.keys():
                print 'error',key,'is missing from mutation data set'
            for mutation in mut_map[protherm[key].pdbwild]:
                if mutation.pkey != key:
                    continue
                try:
                    #print key,protherm[key].pdbwild,mutation.get_pnum(),mutation.chain,protherm[mutation.pkey].ddG,ddg_data[(protherm[key].pdbwild,mutation.chain)][mutation.resnum,mutation.mut]
                    if float(ddg_data[(protherm[key].pdbwild,mutation.chain)][mutation.resnum,mutation.mut]) != []:
                        protherm_ddgs.append(float(protherm[key].ddG))
                        rosetta_ddg = float(ddg_data[(protherm[key].pdbwild,mutation.chain)][mutation.resnum,mutation.mut])
                        rosetta_ddgs.append(rosetta_ddg)

                        stats = STATISTIC_TRACKING(MUTATION_TYPES())
                        stats.add(mutation.wt,mutation.mut)
                        result = ddG_result(key,protherm[key].pdbwild,mutation.resnum,protherm[key].resnum,' ',mutation.chain,mutation.wt,mutation.mut,rosetta_ddg,protherm[key].ddG,stats)
                        results.append(result)
                except:
                    x=None
        #find missing keys
        for key in keys:
            has_key = False
            for result in results:
                if result.key == key:
                    has_key = True
                    break;
            if not has_key:
                print 'missing key',key,'pdbid is',protherm[key].pdbwild
    else:
        benchmark = parse_benchmark_json(args)
        #benchmark = add_rosettanum_to_json(args,benchmark)
        for i in range(0,len(benchmark['data'])-1):
            experiment = benchmark['data'][i]
            stats = STATISTIC_TRACKING(MUTATION_TYPES())
            stats.add_json(experiment['Mutations'][0])
            protherm_ddg = float(experiment['ExperimentalDDGs']['DDG'])
            pdbid = experiment['PDBFileID']
            chain = experiment['Mutations'][0]['Chain']
            ros_num = int(experiment['Mutations'][0]['RosettaNum'])
            mutant_aa = experiment['Mutations'][0]['MutantAA']
            rosetta_ddg = float(ddg_data[(pdbid,chain)][(ros_num,mutant_aa)])
            mutation_num,mutation_icode = split_mutationid(experiment['Mutations'][0]['ResidueID'])
            result = ddG_result(i,pdbid,ros_num,mutation_num,mutation_icode,chain,experiment['Mutations'][0]['WildTypeAA'],mutant_aa,rosetta_ddg,protherm_ddg,stats)
            results.append(result)

    write_mutations(results)
    mutation_types = load_mutation_types()
    
    convert_reu_to_kcal(args,results)
    if args.plottype == 'exact':
        plot_results(args,results,'Prothem ddG','Rosetta ddG',mutation_types)
    if args.plottype == 'cluster' or args.plottype == 'histogram' or args.plottype == 'write outliers':
        original_results = deepcopy(results)
        new_results = cluster_results(results)
        if args.plottype == 'cluster':
            report_error(args,new_results, mutation_types)
        elif args.plottype == 'write outliers':
            write_outliers(args,new_results,2,original_results)

def write_mutations(results):
    with open('targetmutations.txt','w') as tm:
        tm.write('PDBID PDBNUM ROSETTANUM WILD_TYPE MUTATION ROSETTA_DDG PROTHERM_DDG\n')
        for result in results:
            tm.write(result.pdbid+' '+str(result.pdb_num)+' '+str(result.rosetta_num)+' '+result.wt+' '+result.mut+' '+str(result.ros_ddg)+' '+str(result.pro_ddg)+"\n")

def cluster_results(results):
    new_results = []
    for result in results:
        new_result = result
        if result.ros_ddg >= 1:
            new_result.ros_ddg = 1
        if result.ros_ddg <= -1:
            new_result.ros_ddg = -1
        if result.ros_ddg > -1 and result.ros_ddg < 1:
            new_result.ros_ddg = 0
        if result.pro_ddg >= 1:
            new_result.pro_ddg = 1
        if result.pro_ddg <= -1:
            new_result.pro_ddg = -1
        if result.pro_ddg > -1 and result.pro_ddg < 1:
            new_result.pro_ddg = 0
        new_results.append(new_result)
    return new_results

def plot_binding_ddg(args):
    benchmark = parse_benchmark_json(args)
    results = []
    it = 0
    for entry in benchmark['data']:
        mutationnum,mutation_icode = split_mutationid(entry['Mutations'][0]['ResidueID'])
        pdbid = entry['PDBFileID']
        chain = entry['Mutations'][0]['Chain']
        rosnum = entry['Mutations'][0]['RosettaNum']
        wt = entry['Mutations'][0]['WildTypeAA']
        mut = entry['Mutations'][0]['MutantAA']
        proddg = entry['ExperimentalDDGs']['DDG']
        ddgfile = args.jobname+"/"+pdbid+'_'+str(rosnum)+'_'+chain+'_'+mut+"/bindingddg.txt"
        if not os.path.isfile(ddgfile):
            print 'not',ddgfile
            continue
        bindingddg = parse_binding_ddg(ddgfile)/args.scalefactor
        stats = STATISTIC_TRACKING(MUTATION_TYPES())
        stats.add(wt,mut)
        result = ddG_result(it,pdbid,rosnum,mutationnum,mutation_icode,chain,wt,mut,bindingddg,proddg,stats)
        results.append(result)
        it+=1
    mutation_types = load_mutation_types()

    write_mutations(results)
    
    if args.plottype == 'exact':
        plot_results(args,results,'ProTherm ddG','Binding ddG',mutation_types)
    elif args.plottype == 'cluster' or args.plottype == 'write outliers':
        non_clustered = deepcopy(results)
        results = cluster_results(results)
        if args.plottype == 'cluster':
            report_error(args,results,mutation_types)
        elif args.plottype == 'write outliers':
            write_outliers(args,results,2,non_clustered)

#finds the slope of the correlation and uses it to convert reu to kcals
def convert_reu_to_kcal(args,results):
    pro,ros,outlierkeys,firstpearson,oldpearson = filter_outliers(args,results)
    slope = stats.linregress(pro,ros)[0]
    for result in results:
        result.ros_ddg = result.ros_ddg/slope

def load_mutation_types():
    mutation_types = ['small to large','large to small','positive to negative','negative to positive','negative to hydrophobic','hydrophobic to negative','positive to hydrophobic',
            'hydrophobic to positive','non-charged polar to positive','positive to non-charged polar',
            'non-charged polar to negative','negative to non-charged polar','non-charged polar to hydrophobic','hydrophobic to non-charged polar','non-charged polar to non-charged polar',
            'hydrophobic to hydrophobic','charge to charge','involves cysteine',
            'involves proline','same size','everything']
    return mutation_types

def parse_binding_ddg(ddgfile):
    with open(ddgfile,'r') as dfile:
        for line in dfile:
            return float(line.split()[2])

def plot_results(args,results,xaxis,yaxis, mutation_types):

    pp = PdfPages(args.output)
    with open('correlations.txt','w') as cf:
        cf.write('Mutation Type, Pearson\'s R, Pearson\'s R filtered, Removed Entries\n')
        filtered_outliers = []
        for mutation_type in mutation_types:
            in_class_results = []
            for result in results:
                if mutation_type in clasify_mutation(result) or mutation_type == 'everything':
                    in_class_results.append(result)
            protherm_ddgs,rosetta_ddgs,outliers,pearson,newpearson = filter_outliers(args,in_class_results)
            mutations_removed =  len(outliers)
            print mutation_type,str(pearson),str(newpearson),str(mutations_removed)
            cf.write(mutation_type+', '+str(pearson)+', '+str(newpearson)+', '+str(mutations_removed)+'\n')
            #make the plot
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel(xaxis,fontsize=18)
            ax.set_ylabel(yaxis,fontsize=18)
            if args.plottype == 'cluster':
                ax.set_xlim(-1.5,1.5)
                ax.set_ylim(-1.5,1.5)
            ax.scatter(protherm_ddgs,rosetta_ddgs,s=25,c='b')
            ax.set_title(mutation_type,fontsize=24)
            #plot the outliers
            outlier_pro = []
            outlier_ros = []
            filtered_outliers+=outliers
            for outlier in in_class_results:
               if outlier.key not in filtered_outliers:
                   continue
               outlier_pro.append(outlier.pro_ddg)
               outlier_ros.append(outlier.ros_ddg)
            ax.scatter(outlier_pro,outlier_ros,s=25,c='r')
            if args.plottype == 'exact':
                xs,ys = make_plotlines(1.69,0,-10,10,1)
                ax.plot(xs,ys,'r-')
                xs,ys = make_plotlines(1.69,1,-10,10,1)
                ax.plot(xs,ys,'r--')
                xs,ys = make_plotlines(1.69,-1,-10,10,1)
                ax.plot(xs,ys,'r--')
                xs,ys = make_plotlines(1.69,2,-10,10,1)
                ax.plot(xs,ys,'b--')
                xs,ys = make_plotlines(1.69,-2,-10,10,1)
                ax.plot(xs,ys,'b--')
            if args.plottype == 'cluster':
                xs,ys = make_plotlines(1,0,-10,10,1)
                ax.plot(xs,ys,'r-')

            plt.savefig(pp,format='pdf')

    pp.close()

#takes a list of ddG_result objects and returns a list of the rosetta ddg and protherm ddg
def get_ddgs(results):
    ros_ddg = []
    pro_ddg = []
    for result in results:
        ros_ddg.append(result.ros_ddg)
        pro_ddg.append(result.pro_ddg)
    return ros_ddg,pro_ddg

def report_error(args,results,mutation_types):
    pp = PdfPages(args.output)
    with open('errorcount.csv','w') as ec:
        ec.write('mutation class, sameclass, off by 1, off by 2\n')
        for mut_type in mutation_types:
            differences = []
            for result in results:
                if mut_type in clasify_mutation(result) or mut_type == 'everything':
                    differences.append(abs(result.pro_ddg-result.ros_ddg))
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title(mut_type)
            counts = []
            counts.append(differences.count(0))
            counts.append(differences.count(1))
            counts.append(differences.count(2))
            N = 3
            ind = numpy.arange(3)+0.1
            width = 0.35
            ax.bar(ind,counts,width,color='b')
            ax.set_xticks(ind+width-(width/2))
            ax.set_xticklabels(('0','1','2'))
            ax.set_ylim(0,max(counts)+(max(counts)/10))

            plt.savefig(pp,format='pdf')
            ec.write(mut_type+','+str(differences.count(0))+','+str(differences.count(1))+','+str(differences.count(2))+"\n")
    pp.close()

#takes results clasified as 0,1 or 2 and reports the pdbids and numbers for those at the given degree
def write_outliers(args,results,degree,non_clustered):
    with open(args.output,'w') as ol:
        ol.write('PDBID, PDBNUM, ROSETTA NUM, WT, MUT, DEGREE, PRO_DDG, ROS_DDG\n')
        for i in range(0,len(results)-1):
            result = results[i]
            nc_result = non_clustered[i]
            rdegree = abs(result.pro_ddg-result.ros_ddg)
            if rdegree >= degree:
                ol.write(result.pdbid+", "+str(result.pdb_num)+", "+str(result.rosetta_num)+", "+result.wt+", "+result.mut+", "+str(rdegree)+", "+str(nc_result.pro_ddg)+", "+str(nc_result.ros_ddg)+"\n")


#using formula y = mx+b return x and y values
def make_plotlines(m,b,lower,upper,interval):
    x = lower
    ys = []
    xs = []
    while x <= upper:
        y = m*x + b
        ys.append(y)
        xs.append(x)
        x+=interval
    xs = numpy.asarray(xs)
    ys = numpy.asarray(ys)
    return xs,ys

def clasify_mutation(result):
    mutation_classes = []
    if result.stats.sm2l >= 1:
        mutation_classes.append('small to large')
    if result.stats.l2sm >= 1:
        mutation_classes.append('large to small')
    if result.stats.pos2neg >= 1:
        mutation_classes.append('positive to negative')
    if result.stats.neg2pos >= 1:
        mutation_classes.append('negative to positive')
    if result.stats.h2pol >= 1:
        mutation_classes.append('hydrophobic to polar')
    if result.stats.pol2h >= 1:
        mutation_classes.append('polar to hydrophobic')
    if result.stats.pos2ncp >= 1:
        mutation_classes.append('positive to non-charged polar')
    if result.stats.ncp2pos:
        mutation_classes.append('non-charged polar to positive')
    if result.stats.neg2h >=1:
        mutation_classes.append('negative to hydrophobic')
    if result.stats.h2neg >=1:
        mutation_classes.append('hydrophobic to negative')
    if result.stats.pos2h >= 1:
        mutation_classes.append('positive to hydrophobic')
    if result.stats.h2pos >= 1:
        mutation_classes.append('hydrophobic to positive')
    if result.stats.ncp2neg >= 1:
        mutation_classes.append('non-charged polar to negative')
    if result.stats.neg2ncp >= 1:
        mutation_classes.append('negative to non-charged polar')
    if result.stats.ncp2h >= 1:
        mutation_classes.append('non-charged polar to hydrophobic')
    if result.stats.h2ncp >= 1:
        mutation_classes.append('hydrophobic to non-charged polar')
    if result.stats.ncp2ncp >= 1:
        mutation_classes.append('non-charged polar to non-charged polar')
    if result.stats.h2h >= 1:
        mutation_classes.append('hydrophobic to hydrophobic')
    if result.stats.crg2crg >= 1:
        mutation_classes.append('charge to charge')
    if result.stats.cys >= 1:
        mutation_classes.append('involves cysteine')
    if result.stats.pro >= 1:
        mutation_classes.append('involves proline')
    if result.stats.siz2siz >= 1:
        mutation_classes.append('same size')
    return mutation_classes

#swaps the ddg values for the provided publications (PUBMED IDS)
def swap_pmids_protherm(args):
    data = parse_protherm_json(args)
    for mutation in data['data']:
        for experiment in mutation['ExperimentalDDGs']:
            swapped = swap_pmid_ddG(args,experiment)
            if swapped:
                mutation['DDG'] = 'bad'
    print json.dumps(data,sort_keys=True,indent=4,separators=(',',': '))

def swap_pmids_benchmark(args):
    data = parse_benchmark_json(args)
    for mutation in data['data']:
        swapped = swap_pmid_ddG(args,mutation['ExperimentalDDGs'])
        if swapped:
            mutation['DDG'] = 'bad'
    print json.dumps(data,sort_keys=True,indent=4,separators=(',',': '))

#checks the pmid and swaps the ddG if it's in the list for swapped pmids.
def swap_pmid_ddG(args,entry):
    pmid = entry['Publication']
    pmidnum = re.split(':',pmid)[1]
    swapped = False
    if pmidnum in args.pmids:
        entry['DDG'] = -float(entry['DDG'])
        swapped = True
    return swapped


def find_pmids_from_authors(args):
    protherm = parse_protherm(args)
    pmids = []
    for key in protherm.keys():
        entry = protherm[key]
        authors = ''
        reference = ''
        pmid = ''
        for line in entry.lines:
            if line.startswith('AUTHOR'):
                authors = line
            if line.startswith('REFERENCE'):
                if 'PMID:' not in line:
                    continue
                pmid = line.split()[-1]
        for author in args.authors:
            if author.lower().strip() in authors.lower().strip():
                if pmid not in pmids:
                    pmids.append(pmid)
    print pmids

def count_pmids(args):
    print 'counting'
    benchmark = parse_benchmark_json(args)
    pmids = []
    for entry in benchmark['data']:
        pub = entry['ExperimentalDDGs']['Publication']
        pmid = re.split('PMID:',pub)[1]
        if pmid not in pmids:
            pmids.append(pmid)
    for pmid in pmids:
        print pmid
    #print len(pmids)

#one off function to fix the json formatting of the benchmark
def fix_formatting(args):
    benchmark = parse_benchmark_json(args)
    for entry in benchmark['data']:
        try:
            entry['ExperimentalDDGs'] = entry['ExperimentalDDGs'][0]
        except:
            continue
    print json.dumps(benchmark)

#random test function.
def filter_outliers(args,results):
    ros,pro = get_ddgs(results)
    oldpearson = pearsonr(ros,pro)[0]
    firstpearson = deepcopy(oldpearson)
    done = False
    max_filt = 5
    it = 0
    outliers = []
    while not done:
        worst_entry = -1
        worst_index = -1
        for i in range(0,len(results)):
            #don't refilter old results
            if i in outliers:
                continue
            dontneed1,dontneed2,test_pearson = remove_and_pearsons(results,i,outliers)
            if test_pearson > worst_entry:
                worst_entry = test_pearson 
                worst_index = i
        rosddg,proddg,newpearson = remove_and_pearsons(results,worst_index,outliers)
        if newpearson-oldpearson < args.outlier_cutoff:
            done = True
            break
        outliers.append(worst_index)
        it+=1
        if it >= max_filt:
            done = True
        oldpearson = newpearson
    #only changes the stored ddg lists if the pearson'sr changes.
    outlierkeys = []
    if oldpearson != firstpearson:
        ros = rosddg
        pro = proddg
        for outlier in outliers:
            outlierkeys.append(results[outlier].key)
    return pro,ros,outlierkeys,firstpearson,oldpearson

def remove_and_pearsons(results,index,outliers):
    ros,pro = get_ddgs(results)
    toremove = []
    toremove.append(index)
    for outlier in outliers:
        toremove.append(outlier)
    toremove = sorted(toremove,reverse=True)
    for element in toremove:
        del ros[element]
        del pro[element]
    pearsons = pearsonr(ros,pro)[0]
    return ros,pro,pearsons

#starts the binding ddg protocol
def run_binding_ddg(args):
    print 'running binding'
    benchmark = parse_benchmark_json(args)
    if __name__ == '__main__':
        jobs = []
        for i in range(0,len(benchmark['data'])):
            jobs.append((args,i))
        p = multiprocessing.Pool(args.cores,maxtasksperchild=1)
        p.map(que_ddg, jobs)

def que_ddg(inputs):
    args = inputs[0]
    index = inputs[1]
    pdblines = ' '.join(args.pdbs)
    self = os.path.dirname(os.path.realpath(__file__))
    command = ('python %s/ddgstats.py -m "run binding" -bs %s -jn %s -rw %s -ds %s -ss %s -p %s -i %s -r %s\
            '%(self,args.benchmark_set,args.jobname,args.residue_window,args.dumpscript,args.scorescript,pdblines,index,args.relax))
    os.system(command)

def make_singleres_mutfile(entry):
    with open('singlemut.mut','w') as mutfile:
        mutfile.write('total 1\n')
        mutfile.write('1\n')
        wt = entry['Mutations'][0]['WildTypeAA']
        rosnum = entry['Mutations'][0]['RosettaNum']
        mut = entry['Mutations'][0]['MutantAA']
        mutfile.write('%s %s %s'%(wt,rosnum,mut))

def run_separation(args):
    benchmark = parse_benchmark_json(args)
    
    #create directory
    entry = benchmark['data'][args.index]
    mutation = entry['Mutations'][0]
    resnum = entry['Mutations']
    dirname = entry['PDBFileID']+'_'+str(mutation['RosettaNum'])+'_'+mutation['Chain']+'_'+mutation['MutantAA']
    cwd = os.getcwd()
    if not os.path.isdir(args.jobname):
        os.mkdir(args.jobname)
    os.chdir(args.jobname)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)

    # copy files
    pdbfile = ''
    for pdb in args.pdbs:
        if str(entry['PDBFileID']+'_'+mutation['Chain']) in pdb:
            pdbfile = pdb
    shutil.copy(cwd+"/"+pdbfile,'./')
    pdbfile = re.split('/',pdbfile)[-1]
    shutil.copy(cwd+"/"+args.dumpscript,'./')
    shutil.copy(cwd+"/"+args.scorescript,'./')
    if args.relax == 1:
        shutil.copy(cwd+"/runrelax.sh",'./')
        shutil.copy(cwd+"/relax.xml",'./')

    # generate the mutant pdbs
    make_singleres_mutfile(entry)
    command = ('sh %s %s'%(args.dumpscript,pdbfile))
    if len(glob.glob('MUT*')) < 1:
        os.system(command)

    mutpdbs = glob.glob("MUT*")
    if len(mutpdbs) > 1:
        for mpdb in mutpdbs:
            if '.pdb' not in mpdb:
                continue
            if 'extracted' not in mpdb and 'core' not in mpdb:
                mutpdb = mpdb
    else:
        mutpdb = mutpdbs[0]
    if args.relax == 1:
        relax_pdb(mutpdb)
    mutid = pdbtools.get_pdb_id(mutpdb)

    wtpdbs = glob.glob("WT*")
    if len(mutpdbs) > 1:
        for wtpdb in wtpdbs:
            if 'extracted' not in wtpdb and 'core' not in wtpdb:
                wildpdb = wtpdb
    else:
        wildpdb = wtpdbs[0]
    mutid = pdbtools.get_pdb_id(mutpdb)
    
    extract_element(args,wildpdb,mutation['RosettaNum'],'wt')
    extract_element(args,mutpdb,mutation['RosettaNum'],'mut')

    wtpdbid = pdbtools.get_pdb_id(wildpdb)

    allpdbs = [wildpdb,mutpdb,'%s_wt_core.pdb'%wtpdbid,'%s_wt_extracted.pdb'%wtpdbid,'%s_mut_core.pdb'%mutid,'%s_mut_extracted.pdb'%mutid]

    scores = parse_scores('score.sc')
    score_pdbid = ''
    for pdb in allpdbs:
        command = ('sh %s %s'%(args.scorescript,pdb))
        score_pdbid = pdbtools.get_pdb_id(pdb)
        if score_pdbid+'_0001' not in scores.keys():
            os.system(command)
    scores = parse_scores('score.sc')
    print 'calcutlating binding from scores...'
    wt_binding,mut_binding = calculate_binding_ddg(scores)

    with open('bindingddg.txt','w') as bindingfile:
        line = entry['PDBFileID']+' '+str(mutation['ResidueID'])+' '+str(mut_binding-wt_binding)+' '+str(entry['ExperimentalDDGs']['DDG'])
        bindingfile.write(line)


def extract_element(args,pdb,rosnum,tag):
    pdbid = pdbtools.get_pdb_id(pdb)
    residues = pdbtools.get_unopened_residue_list(pdb)
    extracted  = []
    kept = []
    start = rosnum-8
    end = rosnum+8
    for i in range(0,len(residues)):
        residue = residues[i]
        if i < start or i > end:
            kept.append(residue)
        else:
            extracted.append(residue)
    exname = pdbid+"_"+tag+"_extracted.pdb"
    keptname = pdbid+"_"+tag+"_core.pdb"
    pdbtools.write_resis_to_pdb(extracted,exname)
    pdbtools.write_resis_to_pdb(kept,keptname)

def parse_scores(scorefile):
    scores = {}
    if not os.path.isfile(scorefile):
        return scores
    with open(scorefile,'r') as sf:
        for line in sf:
            if line.startswith('SCORE:'):
                line = line.split()
                if line[-1] == 'description':
                    continue
                score = float(line[1])
                pdbid = line[-1]
                scores[pdbid] = score
    return scores

def calculate_binding_ddg(scores):

    for pdbid in scores.keys():
        if 'mut' not in pdbid.lower() and 'wt' not in pdbid:
            fullwt = scores[pdbid]
        if 'MUT' in pdbid and 'core' not in pdbid and 'extracted' not in pdbid:
            fullmut = scores[pdbid]
        if 'wt' in pdbid and 'core' in pdbid:
            wtcore = scores[pdbid]
        if 'wt' in pdbid and 'extracted' in pdbid:
            wtextracted = scores[pdbid]
        if 'mut' in pdbid and 'core' in pdbid:
            mutcore = scores[pdbid]
        if 'mut' in pdbid and 'extracted' in pdbid:
            mutextracted = scores[pdbid]

    #print fullwt,fullmut,wtcore,wtextracted,mutcore,mutextracted

    wt_binding = fullwt - wtextracted
    mut_binding = fullmut - mutextracted
    
    return wt_binding,mut_binding

def collect_binding_ddgs(args):
    benchmark = parse_benchmark_json(args)
    for i in range(0,len(benchmark['data'])):
        entry = benchmark['data'][i]
        mutation = benchmark['data'][0]['Mutations'][0]
        dirname = entry['PDBFileID']+'_'+str(mutation['RosettaNum'])+'_'+mutation['Chain']+'_'+mutation['MutantAA']
        ddgfile = args.jobname+"/"+dirname+"/bindingddg.txt"
        bindingddg = parse_ddg(ddgfile)
        print entry['PDBFileID'],mutation['RosettaNum'],bindingddg,entry['ExperimentalDDGs'][0]['DDG']

def parse_ddg(ddgfile):
    with open(ddgfile,'r') as dfile:
        for line in dfile:
            return line.split()[2]

def get_best_pdb():
    bestpdb = ''
    bestscore = 1e7
    with open('score.sc','r') as scorefile:
        for line in scorefile:
            if line.startswith("SCORE:"):
                line = line.split()
                if line[-1] == 'description':
                    continue
                score = float(line[1])
                if score < bestscore:
                    bestscore = score
                    bestpdb = line[-1]
    return bestpdb

def count_mutations_per_pdb(args):
    benchmark = parse_benchmark_json(args)
    pdbcounts = {}
    for entry in benchmark['data']:
        pdb = entry['PDBFileID']
        if pdb in pdbcounts.keys():
            pdbcounts[pdb]+=1
        else:
            pdbcounts[pdb] = 1

    for pdb in sorted(pdbcounts, key=pdbcounts.get, reverse=False):
        print pdb,pdbcounts[pdb]

def find_two_job_correlation(args):
    benchmark = parse_benchmark_json(args)
    ddg1 = load_ddg_data(args,args.job1)
    ddg2 = load_ddg_data(args,args.job2)
    results1 = []
    results2 = []
    for i in range(0,len(benchmark['data'])-1):
        experiment = benchmark['data'][i]
        stats = STATISTIC_TRACKING(MUTATION_TYPES())
        stats.add_json(experiment['Mutations'][0])
        protherm_ddg = float(experiment['ExperimentalDDGs']['DDG'])
        pdbid = experiment['PDBFileID']
        chain = experiment['Mutations'][0]['Chain']
        ros_num = int(experiment['Mutations'][0]['RosettaNum'])
        mutant_aa = experiment['Mutations'][0]['MutantAA']
        rosetta_ddg1 = float(ddg1[(pdbid,chain)][(ros_num,mutant_aa)])
        rosetta_ddg2 = float(ddg2[(pdbid,chain)][(ros_num,mutant_aa)])
        mutation_num,mutation_icode = split_mutationid(experiment['Mutations'][0]['ResidueID'])
        result1 = ddG_result(i,pdbid,ros_num,mutation_num,mutation_icode,chain,experiment['Mutations'][0]['WildTypeAA'],mutant_aa,rosetta_ddg1,protherm_ddg,stats)
        result2 = ddG_result(i,pdbid,ros_num,mutation_num,mutation_icode,chain,experiment['Mutations'][0]['WildTypeAA'],mutant_aa,rosetta_ddg2,protherm_ddg,stats)
        results1.append(result1)
        results2.append(result2)
    print len(results1),len(results2)


    ddg2exp = list_experiment_values(args,results1)
    ddg1exp = list_experiment_values(args,results2)
    if len(ddg1exp) == 0 or len(ddg2exp) == 0:
        print 'missing ddgs for one or both of the experiments. Exiting'
        exit()
    xys = []
    for i in range(0,len(ddg1exp)):
        xys.append((ddg1exp[i],ddg2exp[2]))
    print len(ddg1exp),len(ddg2exp)
    pearson = pearsonr(ddg1exp,ddg2exp)[0]
    theil = scipy.stats.mstats.theilslopes(xys)
    print pearson,theil


def list_experiment_values(args,results):
    evalues = []
    for result in results:
        if args.plottype not in clasify_mutation(result) and args.plottype != 'everything':
            continue
        experiment = result.ros_ddg
        evalues.append(experiment)
    return evalues

class MUTATION_TYPES:

    def __init__(self):
        self.large = ['F','W','Y','K','R','H','Q','E']
        self.small = ['G','A','V','S','T','C']
        self.polar = ['Y','T','S','H','K','R','E','D','Q','N']
        self.polar_nc = ['Y','T','S','N','Q','H']
        self.positive = ['K','R']
        self.negative = ['E','D']
        self.non_polar = ['F','I','L','V','A','G','M','W']
        self.large_nonpolar = ['F','I','L','V']

class STATISTIC_TRACKING:

    def __init__(self,mut_type):
        self.sm2l = 0
        self.l2sm = 0
        self.pos2neg = 0
        self.neg2pos = 0
        self.h2pol = 0
        self.pol2h = 0
        self.ncp2pos = 0
        self.pos2ncp = 0
        self.ncp2neg = 0
        self.neg2ncp = 0
        self.neg2h = 0
        self.h2neg = 0
        self.pos2h = 0
        self.h2pos = 0
        self.ncp2h = 0
        self.h2ncp = 0
        self.ncp2ncp = 0
        self.h2h = 0
        self.crg2crg = 0
        self.cys = 0
        self.pro = 0
        self.siz2siz = 0
        self.mut_type = mut_type

    #let's you directly use a json entry to add the stats
    def add_json(self,entry):
        wt = entry['WildTypeAA']
        mut = entry['MutantAA']
        self.add(wt,mut)


    def add(self,wt,mut=None):
        #allow this function to take a protherm entry as the first argument
        if  mut is None:
            mut = wt.mutation
            wt = wt.wild

        if wt in self.mut_type.large and mut in self.mut_type.small:
            self.l2sm+=1
        if mut in self.mut_type.large and wt in self.mut_type.small:
            self.sm2l+=1
        if wt in self.mut_type.positive and mut in self.mut_type.negative:
            self.pos2neg+=1
        if wt in self.mut_type.negative and mut in self.mut_type.positive:
            self.neg2pos+=1
        if wt in self.mut_type.non_polar and mut in self.mut_type.polar:
            self.h2pol+=1
        if wt in self.mut_type.polar and mut in self.mut_type.non_polar:
            self.pol2h+=1
        if wt in self.mut_type.polar_nc and mut in self.mut_type.positive:
            self.ncp2pos+=1
        if mut in self.mut_type.polar_nc and wt in self.mut_type.positive:
            self.pos2ncp+=1
        if wt in self.mut_type.polar_nc and mut in self.mut_type.negative:
            self.ncp2neg+=1
        if wt in self.mut_type.negative and mut in self.mut_type.polar_nc:
            self.neg2ncp+=1
        if wt in self.mut_type.negative and mut in self.mut_type.non_polar:
            self.neg2h+=1
        if wt in self.mut_type.non_polar and mut in self.mut_type.negative:
            self.h2neg+=1
        if wt in self.mut_type.positive and mut in self.mut_type.non_polar:
            self.pos2h+=1
        if wt in self.mut_type.non_polar and mut in self.mut_type.positive:
            self.h2pos+=1
        if wt in self.mut_type.polar_nc and mut in self.mut_type.non_polar:
            self.ncp2h+=1
        if wt in self.mut_type.non_polar and mut in self.mut_type.polar_nc:
            self.h2ncp+=1
        if wt in self.mut_type.polar_nc and mut in self.mut_type.polar_nc:
            self.ncp2ncp+=1
        if wt in self.mut_type.non_polar and mut in self.mut_type.non_polar:
            self.h2h+=1
        if (wt in self.mut_type.positive and mut in self.mut_type.positive) or ( wt in self.mut_type.negative and mut in self.mut_type.negative ):
            self.crg2crg+=1
        if (wt in self.mut_type.large and mut in self.mut_type.large) or ( wt in self.mut_type.small and mut in self.mut_type.small ):
            self.siz2siz+=1
        if wt == 'C' or mut == 'C':
            self.cys+=1
        if wt == 'P' or mut == 'P':
            self.pro+=1

    def has_enough(self,amount_wanted,acceptance):
        has_enough = True
        if self.pos2neg < amount_wanted:
            has_enough = False
        if self.neg2pos < amount_wanted:
            has_enough = False
        if self.ncp2pos < amount_wanted:
            has_enough = False
        if self.pos2ncp < amount_wanted:
            has_enough = False
        if self.ncp2neg < amount_wanted:
            has_enough = False
        if self.neg2ncp < amount_wanted:
            has_enough = False
        if self.neg2h < amount_wanted:
            has_enough = False
        if self.h2neg < amount_wanted:
            has_enough = False
        if self.pos2h < amount_wanted:
            has_enough = False
        if self.h2pos < amount_wanted:
            has_enough = False
        if self.ncp2h < amount_wanted:
            has_enough = False
        if self.h2ncp < amount_wanted:
            has_enough = False
        if self.ncp2ncp < amount_wanted:
            has_enough = False
        if self.h2h < amount_wanted:
            has_enough = False
        if self.crg2crg < amount_wanted:
            has_enough = False
        if acceptance == 1:
            return has_enough
        if self.cys < amount_wanted:
            has_enough = False
        if self.pro < amount_wanted:
            has_enough = False
        if acceptance < 3:
            return has_enough
        if self.sm2l < amount_wanted:
            has_enough = False
        if self.l2sm < amount_wanted:
            has_enough = False
        if self.h2pol < amount_wanted:
            has_enough = False
        if self.pol2h < amount_wanted:
            has_enough = False
        if self.siz2siz < amount_wanted:
            has_enough = False

        return has_enough

    def added_new(self,stats,amount_wanted,acceptance=3):
        add_new = False
        if self.pos2neg > stats.pos2neg and self.pos2neg <= amount_wanted:
            add_new = True
        if self.neg2pos > stats.neg2pos and self.neg2pos <= amount_wanted:
            add_new = True
        if self.ncp2pos > stats.ncp2pos and self.ncp2pos <= amount_wanted:
            add_new = True
        if self.pos2ncp > stats.pos2ncp and self.pos2ncp <= amount_wanted:
            add_new = True
        if self.neg2h > stats.neg2h and self.neg2h <= amount_wanted:
            add_new = True
        if self.h2neg > stats.h2neg and self.h2neg <= amount_wanted:
            add_new = True
        if self.pos2h > stats.pos2h and self.pos2h <= amount_wanted:
            add_new = True
        if self.h2pos > stats.h2pos and self.h2pos <= amount_wanted:
            add_new = True
        if self.ncp2neg > stats.ncp2neg and self.ncp2neg <= amount_wanted:
            add_new = True
        if self.neg2ncp > stats.neg2ncp and self.neg2ncp <= amount_wanted:
            add_new = True
        if self.ncp2h > stats.ncp2h and self.ncp2h <= amount_wanted:
            add_new = True
        if self.h2ncp > stats.h2ncp and self.h2ncp <= amount_wanted:
            add_new = True
        if self.ncp2ncp > stats.ncp2ncp and self.ncp2ncp <= amount_wanted:
            add_new = True
        if self.h2h > stats.h2h and self.h2h <= amount_wanted:
            add_new = True
        if self.crg2crg > stats.crg2crg and self.crg2crg <= amount_wanted:
            add_new = True
        if acceptance == 1:
            return add_new
        if self.cys > stats.cys and self.cys <= amount_wanted:
            add_new = True
        if self.pro > stats.pro and self.pro <= amount_wanted:
            add_new = True
        if acceptance < 3:
            return add_new
        if self.sm2l > stats.sm2l and self.sm2l <= amount_wanted:
            add_new = True
        if self.l2sm > stats.l2sm and self.l2sm <= amount_wanted:
            add_new = True
        if self.h2pol > stats.h2pol and self.h2pol <= amount_wanted:
            add_new = True
        if self.pol2h > stats.pol2h and self.pol2h <= amount_wanted:
            add_new = True
        if self.siz2siz > stats.siz2siz and self.siz2siz <= amount_wanted:
            add_new = True

        return add_new

    def to_dict(self):
        stats = {}
        stats['small to large'] = self.sm2l
        stats['large to small'] = self.l2sm
        stats['positive to negative'] = self.pos2neg
        stats['negative to positive'] = self.neg2pos
        stats['hydrophobic to polar'] = self.h2pol
        stats['polar to hydrophobic'] = self.h2pol
        stats['non-charged polar to positive'] = self.ncp2pos
        stats['positive to non-charged polar'] = self.pos2ncp
        stats['non-charged polar to negative'] = self.ncp2neg
        stats['negative to non-charged polar'] = self.neg2ncp
        stats['negative to hydrophobic'] = self.neg2h
        stats['hydrophibic to negative'] = self.h2neg
        stats['positive to hydrophobic'] = self.pos2h
        stats['hydrophobic to postiive'] = self.h2pos
        stats['non-charged polar to hydrophobic'] = self.ncp2h
        stats['hydrophobic to non-charged polar'] = self.h2ncp
        stats['non-charged polar to non-charged polar'] = self.ncp2ncp
        stats['hydrophobic to hydrophobic'] = self.h2h
        stats['like charge to like charge'] = self.crg2crg
        stats['involves cysteine'] = self.cys
        stats['involves proline'] = self.pro
        stats['same size to same size'] = self.siz2siz
        return stats

class PROTHERM_ENTRY:

    def __init__(self,entrynum,all_lines):
        self.entrynum = entrynum
        self.lines = all_lines
        self.wild = None
        self.mutation = None
        self.resnum = ''

    def parse_lines(self):
        self.ddG = ''
        for line in self.lines:
            if line.startswith('PDB_wild'):
                self.pdbwild = get_line_data(line)
            if line.startswith('PDB_mutant'):
                self.pdbmut = get_line_data(line)
            if line.startswith('ddG_H2O'):
                lline = line.split()
                if len(lline) > 1:
                    try:
                        self.ddG = -float(lline[1])
                    except:
                        self.ddG = ''
            if line.startswith('ddG') and self.ddG == '':
                lline = line.split()
                if len(lline) > 1:
                    try:
                        self.ddG = -float(lline[1])
                    except:
                        self.ddG = ''
                else:
                    self.ddG = ''

            if line.startswith('MUTATION'):
                lline = line.split()
                if len(lline)> 1:
                    wild = lline[1]
                    if wild != 'wild':
                        self.wild = wild
                        if(len(lline) > 4):
                            self.mutation = 'multi'
                        elif len(lline) == 4:
                            self.mutation = lline[3]
                            self.resnum = lline[2]
            if line.startswith("MUTATED_CHAIN"):
                chain = get_line_data(line)
                if chain == '-' or chain == '_':
                    chain = 'HD'
                #if chain != 'HD' and chain != 'A' and chain.strip() != '':
                #    print chain
                self.mutchain = chain

class Mutation:

    def __init__(self,resnum,chain,wt,mut):
        self.resnum = resnum
        self.chain = chain
        self.wt = wt
        self.mut = mut

    #gives the Mutation object the protherm key it came from
    def set_protherm_key(self,key):
        self.pkey = key

    #returns the protherm key
    def get_protherm_key(self):
        try:
            return self.pkey
        except:
            print 'No key associated with this mutation'
            return None
    
    #getters and setters for the number in pdb numbering format (as opposed to Rosetta numbering)
    def set_pnum(self,pnum):
        self.pnum = pnum

    def get_pnum(self):
        try:
            return self.pnum
        except:
            print 'No pnumber has been set for this mutation'
            return None


class ddG_result:

    def __init__(self,key,pdbid,rosetta_num,pdb_num,pdb_icode,chain,wt,mut,ros_ddg,pro_ddg,stats):
        self.key = key
        self.pdbid = pdbid
        self.rosetta_num = rosetta_num
        self.pdb_num = pdb_num
        self.chain = chain
        self.wt = wt
        self.mut = mut
        self.ros_ddg = ros_ddg
        self.pro_ddg = pro_ddg
        self.stats = stats

class outlier:

    def __init__(self,pdbid,pdbnum,rosnum,wt,mut,rosddg,proddg):
        self.pdbid = pdbid
        self.pdbnum = pdbnum
        self.rosnum = rosnum
        self.wt = wt
        self.mut = mut
        self.rosddg = rosddg
        self.proddg = proddg

    def print_outlier(self):
        print self.pdbid,self.pdbnum,self.rosnum,self.wt,self.mut,self.rosddg,self.proddg
