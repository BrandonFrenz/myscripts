#!/usr/bin/python
import argparse
import numpy as np
from collections import defaultdict
import random
import pdbtools
import math
import re

def main():
    args = parseargs()
    args = parse_resis(args)
    run_apply(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--coordinates',help='The coordinate file')
    parser.add_argument('-dc','--distancecut',default=4.3,type=float,help='The distance cutoff for all points')
    parser.add_argument('-o','--outpdb',default='pdbfrompoints.pdb',help='The name of the output pdb')
    parser.add_argument('-lr','--lower_res',help='The lower residue')
    parser.add_argument('-ur','--upper_res',help='The upper residue')
    parser.add_argument('-p','--pdb',help='The pdb to use for lower and upper')
    args = parser.parse_args()
    return args

def parse_resis(args):
    num1 = int(re.split('\D+',args.lower_res)[0])
    chain1 = re.split('\d+',args.lower_res)[1]
    respair1 = (num1,chain1)
    num2 = int(re.split('\D+',args.upper_res)[0])
    chain2 = re.split('\d+',args.upper_res)[1]
    respair2 = (num2,chain2)
    args.lower_res = respair1
    args.upper_res = respair2
    return args

def run_apply(args):

    resis = pdbtools.get_unopened_residue_list(args.pdb)
    res1 = get_resi(resis,args.lower_res)
    res2 = get_resi(resis,args.upper_res)

    point1 = residue_to_point(res1)
    point2 = residue_to_point(res2)

     
    mindist = res2.num-res1.num*2.7
    maxdist = res2.num-res1.num*4.3
    print maxdist
    exit()
    
    coords = get_coordinates(args)
    coords = clean_coordinates(args,coords)
    neighbor_graph = find_neighbors(args,coords)
    #verify_neighbors(neighbor_graph)
    #exit()
    
    cpstart = find_closest_point(point1,neighbor_graph)
    cpend = find_closest_point(point2,neighbor_graph)
    print cpstart,cpend
    total_dist = random.uniform(mindist,maxdist)
    lower_points,lower_dist,upper_points,upper_dist = initialize_loops(cpstart,cpend,neighbor_graph,total_dist)
    for i in range(1,250):
        lower,upper,dist = run_montecarlo(neighbor_graph,lower_points,lower_dist,upper_points,upper_dist,mindist,maxdist)
        print dist
        if dist < 6:
            name = 'points_'+str(dist)+'_'+str(i)+'.pdb'
            pdb_from_points(name,neighbor_graph,lower,upper)


def run_montecarlo(neighbor_graph,lower_points,lower_dist,upper_points,upper_dist,mindist,maxdist):
    
    closest_gap = gap_distance(lower_points,upper_points,neighbor_graph)
    current_gap = closest_gap
    best_lower = lower_points
    best_upper = upper_points
    
    montecount = 0
    kt = 200
    for it in range(0,7):
        kt = kt/2
        for i in range(0,100):
            total_dist = random.uniform(mindist,maxdist)
            points_to_sample = []
            lower = random.choice(['True','False'])
            newdist = 0
            startpoint = 0
            if lower:
                points_to_sample = lower_points
                newdist = total_dist-upper_dist
                startpoint = random.randint(1,max(1,len(lower_points)-2))
            else:
                points_to_sample = upper_points
                newdist = total_dist-lower_dist
                startpoint = random.randint(1,max(1,len(upper_points)-2))
        
            #sample new points
            newpoints,newdist = resample_points(points_to_sample,newdist,startpoint,neighbor_graph)
            if newdist < closest_gap:
                closest_gap = newdist
                best_lower = lower_points
                best_upper = upper_points
            
            accept = False
            arg = (newdist-current_gap)/kt
            probability = math.exp(-arg)
            if probability > 1.0:
                probability = 1
            if probability >= random.uniform(0,1):
                accept = True
            if accept:
                current_gap = newdist
                if lower:
                    lower_points = newpoints
                    lower_dist = newdist
                else:
                    upper_points = newpoints
                    upper_dist = newdist
    
    return best_lower,best_upper,closest_gap
    
def get_coordinates(args):
    coordinates = []
    with open(args.coordinates) as coordfile:
        for line in coordfile:
            if line.startswith('ATOM'):
                line = line.split()
                x = float(line[2])
                y = float(line[3])
                z = float(line[4])
                coordinates.append((x,y,z))
    return coordinates

def clean_coordinates(args,coords):
    cleancoords = []
    for coord in coords:
        hasneighbor = False
        for coord2 in coords:
            if coord == coord2:
                continue
            dist = point_distance(coord,coord2)
            if dist < args.distancecut:
                hasneighbor = True
        if hasneighbor == True:
            cleancoords.append(coord)
        else:
            print coord,'has no neighbors'
    return cleancoords

def find_neighbors(args,coordinates):
    neighborgraph = defaultdict(dict)
    it = 1
    for coord in coordinates:
        neighbors = defaultdict(list)
        nit = 1
        mindist = 0
        for ncoord in coordinates:
            if coord == ncoord:
                continue
            dist = point_distance(coord,ncoord)
            if dist < mindist or mindist == 0:
                mindist = dist
            if dist < args.distancecut:
                neighbors[nit] = dist
            nit+=1
        #skip any coordinates that have no neighbors they are a dead end
        if len(neighbors.keys()) == 0:
            print 'skipping coord',it,'at',coord,'for not having any neighbors'
            ##it+=1
            continue
        #print it,len(neighbors)
        neighborgraph[it]['neighbors'] = neighbors
        neighborgraph[it]['coordinate'] = coord
        it+=1
    return neighborgraph

def point_distance(point1,point2):
    fc = np.array(point1)
    sc = np.array(point2)
    dist = np.linalg.norm(fc-sc)
    return dist

def gap_distance(lower_points,upper_points,neighbor_graph):
    dist = point_distance(neighbor_graph[lower_points[-1]]['coordinate'],neighbor_graph[upper_points[-1]]['coordinate'])
    return dist

def find_closest_point(coord,neighbormap):
    closestpoint = ''
    closestdist = 0
    for point in neighbormap.keys():
        dist = point_distance(coord,neighbormap[point]['coordinate'])
        if dist < closestdist or closestdist == 0:
            closestdist = dist
            closestpoint = point
    return closestpoint

def initialize_loops(start,end,neighborgraph,total_dist):
    lower = random.uniform(0,total_dist)
    upper = total_dist-lower
    lowerpoints = [start]
    upperpoints = [end]
    lower_points,ldist = resample_points(lowerpoints,lower,0,neighborgraph)
    upper_points,udist = resample_points(upperpoints,upper,0,neighborgraph)
    return lower_points,ldist,upper_points,udist

def resample_points(points,distance,start,neighbor_graph):
    newpoints = []
    for i in range(0,min(len(points),start+1)):
        print len(points)
        print i
        newpoints.append(points[i])
    currentdist = 0
    for i in range(0,len(newpoints)-1):
        currentdist+=neighbor_graph[newpoints[i]]['neighbors'][newpoints[i+1]]
    if currentdist > distance:
        return newpoints,currentdist
    while currentdist < distance:
        rint = random.randint(0,len(neighbor_graph[newpoints[-1]]['neighbors'].keys())-1)
        cp = newpoints[-1]
        np = sorted(neighbor_graph[cp]['neighbors'].keys())[rint]
        newpoints.append(np)
        currentdist+=neighbor_graph[newpoints[-2]]['neighbors'][np]
    return newpoints,currentdist

#def pdb_from_multi_points(name,neighbor_graph,lower_points,upper_points):
#    combined_points = lower_points+upper_points
#    pdb_from_points(name,neighbor_graph,combined_points)

def pdb_from_points(name,neighbor_graph,lower,upper):
    resis = []
    for point in lower:
        resi = point_to_residue(neighbor_graph[point]['coordinate'],point,'A')
        resis.append(resi)
    for point in upper:
        resi = point_to_residue(neighbor_graph[point]['coordinate'],point,'B')
        resis.append(resi)
    pdbtools.write_resis_to_pdb(resis,name)

def point_to_residue(point,resnum,chain):
    atom = pdbtools.Atom('ATOM  ',resnum,' MG ',' ',' ',point[0],point[1],point[2],' ',' ',' ','Mg',' ')
    resi = pdbtools.Residue(resnum,chain,'UNK',[atom])
    resi.isterm = False
    return resi

def residue_to_point(resi):
    ca = pdbtools.get_ca(resi)
    point = (ca.x,ca.y,ca.z)
    return point

def verify_neighbors(neighbor_graph):
    for key in sorted(neighbor_graph.keys()):
        #print key
        for neighbor in neighbor_graph[key]['neighbors']:
            test = True
        test = neighbor_graph[key]['coordinate']

def get_resi(residues,num_chain):
    for res in residues:
        if res.num == num_chain[0]:
            if res.chain == num_chain[1]:
                return res
    print 'no res for',num_chain,'exiting'
    exit()

main()
