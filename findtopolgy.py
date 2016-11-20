#!/usr/bin/python
import argparse
import numpy as np
from collections import defaultdict
import random
import pdbtools
import math

def main():
    args = parseargs()
    run_apply(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--coordinates',help='The coordinate file')
    parser.add_argument('-dc','--distancecut',default=4.3,type=float,help='The distance cutoff for all points')
    parser.add_argument('-o','--outpdb',default='pdbfrompoints.pdb',help='The name of the output pdb')
    args = parser.parse_args()
    return args

def run_apply(args):
    
    #these magic numbers are just for testing
    point = (296,-61,141)
    point2 = (267,-66,149)
    mindist = 55*2.7
    maxdist = 55*4.3
    
    coords = get_coordinates(args)
    neighbor_graph = find_neighbors(args,coords)
    
    cpstart = find_closest_point(point,neighbor_graph)
    cpend = find_closest_point(point2,neighbor_graph)
    total_dist = random.uniform(mindist,maxdist)
    lower_points,lower_dist,upper_points,upper_dist = initialize_loops(cpstart,cpend,neighbor_graph,total_dist)
    for i in range(1,250):
        lower,upper,dist = run_montecarlo(neighbor_graph,lower_points,lower_dist,upper_points,upper_dist,mindist,maxdist)
        if dist < 6:
            name = 'points_'+str(dist)+'_'+str(i)+'.pdb'
            pdb_from_multi_points(name,neighbor_graph,lower,upper)


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
                startpoint = random.randint(1,max(1,len(lower_points)-1))
            else:
                points_to_sample = upper_points
                newdist = total_dist-lower_dist
                startpoint = random.randint(1,max(1,len(upper_points)-1))
        
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
            print 'skipping coord',nit,'at',coord,'for not having any neighbors'
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
    for i in range(0,start+1):
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

def pdb_from_multi_points(name,neighbor_graph,lower_points,upper_points):
    combined_points = lower_points+upper_points
    pdb_from_points(name,neighbor_graph,combined_points)

def pdb_from_points(name,neighbor_graph,points):
    resis = []
    for point in points:
        resi = point_to_residue(neighbor_graph[point]['coordinate'],point)
        resis.append(resi)
    pdbtools.write_resis_to_pdb(resis,name)

def point_to_residue(point,resnum):
    atom = pdbtools.Atom('ATOM  ',resnum,' MG ',' ',' ',point[0],point[1],point[2],' ',' ',' ','Mg',' ')
    resi = pdbtools.Residue(resnum,'A','UNK',[atom])
    resi.isterm = False
    return resi

main()
