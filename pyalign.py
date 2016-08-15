#!/usr/bin/python

import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
import sys, time, os
import pymol
import argparse
 

#Aligns the mobile and static structure using pymol align
def alignstructures(mobile,static):
    pymol.finish_launching()
    staticStructurePath = os.path.abspath(static)
    staticStructureName = staticStructurePath.split('/')[-1].split('.pdb')[0]
    mobileStructurePath = os.path.abspath(mobile)
    mobileStructureName = mobileStructurePath.split('/')[-1].split('.pdb')[0]
 
    # Load Structures
    pymol.cmd.load(mobileStructurePath, mobileStructureName)
    pymol.cmd.load(staticStructurePath, staticStructureName)
 
    pymol.cmd.align(mobileStructureName,staticStructureName,1)
    #time.sleep(1) # Dunno why, but if I don't wait, structures do not align properly..
    # Save Superimposition
    # save(file, selection, state (0 default), format)
    #pymol.cmd.save("new%s.pdb" %(mobileStructureName), mobileStructureName, 0, 'pdb')
    pymol.cmd.save("a_%s.pdb" %(mobileStructureName), mobileStructureName, 0, 'pdb')

    pymol.cmd.quit()
 
