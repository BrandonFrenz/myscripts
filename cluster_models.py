#!/usr/bin/python

import os,sys
import fileinput
from glob import glob
from math import sqrt
from numpy import average, std
import argparse
import pdbtools

#pdblist = glob('compoundB_interface_*.pdb')
#pdblist.sort()

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structures', nargs="+", help='The structures to cluster')
    parser.add_argument('-ct', '--cluster_tolerance', type=float, default=1.0, help='The cluster tolerance for the ligands')
    args = parser.parse_args()
    return args

def get_atom(filename, resname, atomname):
  for line in fileinput.input(filename):
    if line.startswith('HETATM') or line.startswith('ATOM'):
      line_res = line[17:20].strip()
      line_atom = line[12:16].strip()
      if (line_res == resname) and (line_atom == atomname):
        xval = float(line[31:38]); yval = float(line[38:46]); zval = float(line[46:54])
        coords = [xval,yval,zval]
  return(coords)


def get_relevant_atoms(filename, chain_name):
  backbone = []
  atom_list = ['C','CA','N']
  for line in fileinput.input(filename):
    if line.startswith('ATOM'):
      if (line[21] == chain_name) and (line[12:16].strip() in atom_list):
        xval = float(line[31:38]); yval = float(line[38:46]); zval = float(line[46:54])
        coords = [xval,yval,zval]
        backbone.append(coords)
  return(backbone)


def rmsd(coords1,coords2):
  if len(coords1) != len(coords2):
    sys.exit('error: cannot compute rmsd between coordinate sets of unequal length')
  dist_sum = 0
  for x in range(len(coords1)):
    dist_sum += xyzdist(coords1[x],coords2[x])
  return( sqrt( dist_sum/float(len(coords1)) ) )


def xyzdist(p1,p2):
  dist = sqrt( (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2 )
  return dist

def parse_scores(pdbfiles):
    scores = {}
    for pdbfile in pdbfiles:
        score = float(pdbfile.split('_')[1].split('.pdb')[0])
        scores[pdbfile] = score
    return scores

def sort_pdbs(structure_list):
    pairs = []
    for pdb in structure_list:
        modnum = int(pdb.split('rank')[1].split('_')[0])
        pairs.append((modnum, pdb))
    return [x[1] for x in sorted(pairs)]


def main():

    args = parseargs()
    scores = parse_scores(args.structures)
    pdblist = sort_pdbs(args.structures)

    # first loop through pdbs to gather coordinates
    all_coords = []
    for pdb in pdblist:
      modnum = int(pdb.split('rank')[1].split('_')[0])
      modscore = scores[pdb]
      ligandatoms = []
      ligand = pdbtools.get_unopened_residue_list(pdb)[-1]
      ligandcoords = []
      for atom in ligand.atoms:
          if atom.element != "H":
              ligandcoords.append([atom.x, atom.y, atom.z])
      all_coords.append([ligandcoords,modnum,modscore,pdb])
    
    clusters = []
    ### perform agglomerative clustering
    for x in range(0,len(all_coords)):
      ### add first element as first cluster
      added = 0
      if len(clusters) == 0:
        clusters.append([all_coords[x]])
        added = 1
      else:
        ### compute rmsd of coords to all existing clusters
        ### if rmsd within cluster_tol, add to existing cluster
        for cluster in clusters:
          for mod in cluster:
    #      for mod in cluster[:1]:    ### if you only want to use the head of the cluster
            clust_rmsd = rmsd( mod[0], all_coords[x][0] )
    #        print('rmsd between %s and %s = %.2f'%(mod[1], all_coords[x][1], clust_rmsd))
            if clust_rmsd <= args.cluster_tolerance:
              cluster.append( all_coords[x] )
              added = 1
              break
      ### if current model doesn't fit in existing cluster, create a new one
      if not added:
        clusters.append( [all_coords[x]] )
    
    fout = open('cluster_data.txt','w')
    cluster_string = ''
    for x, cluster in enumerate(clusters):
      cluster_string = ' '
      cluster_scores = []
      max_linker_distances = []
      min_linker_distances = []
      for y, mod in enumerate(cluster):
        cluster_scores.append(mod[2])
        cluster_string = cluster_string+str(mod[1])+' '
        if (y+1)%10 == 0:
          cluster_string = cluster_string[:-1]+'\n'+'               '
      line1 = 'cluster %2i has %3i members:'%(x+1,len(cluster))
      line2 = '    average score = %.2f ( %.2f )'%(average(cluster_scores),std(cluster_scores))
      line5 = '    model(s): %s\n'%(cluster_string[:-1])
      print(line1)
      print(line2)
      print(line5)
    #  print(cluster_string[:-1])
      fout.write(line1+'\n')
      fout.write(line2+'\n')
      fout.write(line5+'\n')
    #  fout.write(cluster_string[:-1]+'\n')
    
    #  print('rmsd between models %3i and %3i = %.2f'%(1,x+1,rmsd(all_coords[0][0],all_coords[x][0])))

if __name__ == '__main__':
    main()
