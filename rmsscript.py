#!/usr/bin/python
import argparse

def main():
    args = parseargs()
    get_rms(args)

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p1','--pdb1',help='The first pdb')
    parser.add_argument('-n','--native',help='The native pdb')
    args = parser.parse_args()
    return args

def parse_file(filename):
    points = {}
    with open(filename,'w') as coordfile:
        for line in coordfile:
            data = line.split()
            resnum = int(data[0])
            x = float(data[1])
            y = float(data[2])
            z = float(data[3])
            points[resnum] = (x,y,z)
    return points

def get_rms(args):
    consensus = parse_file(args.pdb1)
    native = parse_file(args.native)

    #consensus_list = []
    #native_list = []
    total_dst_squared = 0
    for resnum in consensus.keys():
        #consensus_list.append(consensus[resnum])
        #native_list.append(native[resnum])
        consensus_coord = np.array(consensus[resnum])
        native_coord = np.array(native[resnum])
        dist = np.linalg.norm(firstcoords-secondcoords)
        total_dst_squared += dist*dist
    rms = np.sqrt(total_dst_squared/len(consensus.keys())



main()
