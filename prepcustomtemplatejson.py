#!/usr/bin/env python
import argparse
import json
import os

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--structures', nargs='+', help='The structure files')
    args = parser.parse_args()
    return args

def main():
    args = parseargs()
    mylist = []
    count=1
    for structure in args.structures:
        base = os.path.basename(structure)
        pdbid = os.path.splitext(base)[0]
        template_dict = {}
        template_dict['allow-ligands'] = "false"
        template_dict['pdb-key'] = f'brandon/for-karen/hm-with-dna/{base}'
        template_dict['weight'] = 10
        template_dict['templateID'] = pdbid
        mylist.append(template_dict)
        count+=1
    print(json.dumps(mylist))


if __name__ == '__main__':
    main()
