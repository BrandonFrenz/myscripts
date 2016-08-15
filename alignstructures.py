#!/usr/bin/python


import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
import sys, time, os
import pymol
import argparse
 
pymol.finish_launching()

def main():
    
    args = parseargs()
    alignstructures(args)
 

#read the command line arguments 
def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--static', help="the static, non-moving, structure")
    parser.add_argument('-m', '--mobiles', nargs="+", help="the list of moving structures")
    args = parser.parse_args()
    return args

#Align all of the mobile structures to the static and save
def alignstructures(args):
    staticStructurePath = os.path.abspath(args.static)
    staticStructureName = staticStructurePath.split('/')[-1].split('.pdb')[0]
    for mobilestruc in args.mobiles:
        mobileStructurePath = os.path.abspath(mobilestruc)
        mobileStructureName = mobileStructurePath.split('/')[-1].split('.pdb')[0]
 
        # Load Structures
        pymol.cmd.load(mobileStructurePath, mobileStructureName)
        pymol.cmd.load(staticStructurePath, staticStructureName)
 
        pymol.cmd.align(mobileStructureName,staticStructureName,1)
        #time.sleep(1) # Dunno why, but if I don't wait, structures do not align properly..
        # Save Superimposition
        # save(file, selection, state (0 default), format)
        #pymol.cmd.save("new%s.pdb" %(mobileStructureName), mobileStructureName, 0, 'pdb')
        pymol.cmd.save(mobileStructurePath, mobileStructureName, 0, 'pdb')


#def alignligands(mobile,static):
#    # Load Structures
#    pymol.cmd.load(mobileStructurePath, mobileStructureName)
#    pymol.cmd.load(staticStructurePath, staticStructureName)

main()

pymol.cmd.quit()
 
