#!/usr/bin
import re
import amino_acids
import os
import numpy as np

def get_pdb_id(pdbpath):
    pdbid = re.split('/|.pdb',pdbpath)[-2]
    return pdbid

def convert_resis_to_ala(pdb,residuestochange):
    newpdb = []
    backbones = ['N','CA','C','O','CB']
    for line in pdb:
        if line.startswith('ATOM'):
            ll = list(line)
            restype = ''.join(ll[17:20])
            resnum = int(''.join(ll[23:26]).strip())
            if resnum not in residuestochange:
                newpdb.append(line)
                continue
            if restype in amino_acids.longer_names and restype != 'GLY' and restype != 'ALA':
                atomtype = "".join(ll[12:16]).strip()
                if atomtype not in backbones:
                    continue
                else:
                    line = re.sub(restype,'ALA',line)
                    newpdb.append(line)
            else:
                newpdb.append(line)
        else:
            newpdb.append(line)
    
    return newpdb

def convert_resis(pdb,residuestochange,newrestype):
    newpdb = []
    backbones = ['N','CA','C','O','CB']
    for line in pdb:
        if line.startswith('ATOM'):
            ll = list(line)
            restype = ''.join(ll[17:20])
            resnum = int(''.join(ll[23:26]).strip())
            if resnum not in residuestochange:
                newpdb.append(line)
                continue
            if restype in amino_acids.longer_names:
                atomtype = "".join(ll[12:16]).strip()
                if atomtype not in backbones:
                    continue
                else:
                    line = re.sub(restype,newrestype,line)
                    newpdb.append(line)
            else:
                newpdb.append(line)
        else:
            newpdb.append(line)
    
    return newpdb

def get_seq_from_resis(residues):
    seq = []
    for residue in residues:
        if residue.name not in amino_acids.longer_names:
            seq.append('X')
        else:
            seq.append(amino_acids.longer_names[residue.name])
    return ''.join(seq)

def write_pdb(pdbfile,name):
    with open(name,'w') as newpdb:
        for line in pdbfile:
            newpdb.write(line)

def atom_from_pdbline(line):
    line = line.strip()
    ll = list(line)
    record = ''.join(ll[0:6])
    num = int(''.join(ll[6:11]))
    atomid = ''.join(ll[12:16])
    ali = ll[16]
    achar = ll[26]
    x = float(''.join(ll[30:38]))
    y = float(''.join(ll[39:46]))
    z = float(''.join(ll[47:54]))
    occupancy = ''.join(ll[55:60])
    temp = ''.join(ll[60:66])
    segid = ''.join(ll[72:76])
    element = ''.join(ll[76:78])
    charge = ''.join(ll[79:80])
    atom = Atom(record,num,atomid,ali,achar,x,y,z,occupancy,temp,segid,element,charge)
    return atom

def get_residue_list(pdbfile):
    residues = []
    previousresid = ('','')
    currentres = ''
    linecount = 0
    for line in pdbfile:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            ll = list(line)
            resnum = int(''.join(ll[22:26]))
            chainid = ''.join(ll[21])
            currentresid = (resnum,chainid)
            if previousresid != currentresid:   
                previousresid = currentresid
                if currentres != '':
                    isterm = False
                    if pdbfile[linecount-1].startswith('TER'):
                        isterm = True
                    currentres.isterm = isterm
                    residues.append(currentres)
                resname = ''.join(ll[17:20])
                atom = atom_from_pdbline(line)
                currentres = Residue(resnum,chainid,resname,[])
                currentres.atoms.append(atom)
            else:
                atom = atom_from_pdbline(line)
                currentres.atoms.append(atom)
        linecount+=1
    currentres.isterm = True
    residues.append(currentres)
    return residues

def count_cas(residues):
    cacount = 0
    for resi in residues:
        for atom in resi.atoms:
            if atom.atomid == ' CA ':
                cacount+=1
    return cacount

def get_center_of_mass(residues):
    total_cas = count_cas(residues)
    cax = 0
    cay = 0
    caz = 0
    for residue in residues:
        ca = residue.ca()
        cax += ca.x/total_cas
        cay += ca.y/total_cas
        caz += ca.z/total_cas
    return [cax,cay,caz]


#this function will result in ter statements being added whereever non continuous residue numbering occurs
def add_ters_to_noncontres(residues):
    updated_residues = []
    for it in range(0, len(residues)-1):
        res = residues[it]
        if res.num != residues[it+1].num-1:
            res.isterm = True
        updated_residues.append(res)
    return updated_residues

def make_pdblines_from_residues(reslist):
    pdbfile = []
    previousres = (0,'0')
    reslist = sorted(reslist, key=lambda x: (x.chain,x.num))
    for res in reslist:
        previousres = (res.num,res.chain)
        for atom in res.atoms:
            line = list(' '*80)
            line[0:6] = atom.record
            line[6:11] = ' '*(5-len(str(atom.num)))+str(atom.num)
            line[12:16] = atom.atomid
            line[16] = atom.ali
            line[17:20] = res.name
            line[21] = res.chain
            line[22:26] = ' '*(4-len(str(res.num)))+str(res.num)
            line[26] = atom.Acode
            line[30:38] = format_pdb_coord(atom.x)
            line[38:46] = format_pdb_coord(atom.y)
            line[46:54] = format_pdb_coord(atom.z)
            line[55:60] = atom.occupancy
            line[60:66] = atom.tempfact
            line[72:76] = atom.segid
            line[76:78] = atom.element
            line[78:80] = atom.charge
            line = ''.join(line)+'\n'
            pdbfile.append(line)
        #add a ter statement if the residue numbering increased by more than 1
        if res.isterm:
            terline = 'TER\n'
            pdbfile.append(terline)
    return pdbfile

def format_pdb_coord(coord):
    coord = str(coord).split('.')
    predec = coord[0]
    postdec = coord[1]
    newstr = ' '*(4-len(predec))+predec+'.'+postdec+'0'*(3-len(postdec))
    return newstr

def atomlist_rms(atoms1,atoms2):
    assert(len(atoms1) == len(atoms2)), 'unequal number of atoms exiting'
    n = len(atoms1)
    firstcoords = []
    for atom in atoms1:
        firstcoords.append([atom.x,atom.y,atom.z])

    secondcoords = []
    for atom in atoms2:
        secondcoords.append([atom.x,atom.y,atom.z])
    firstcoords = np.array(firstcoords)
    secondcoords = np.array(secondcoords)
    rms = np.linalg.norm(firstcoords-secondcoords)/np.sqrt(n)
    return rms

def atom_dist(atom1,atom2):
    firstcoords = []
    firstcoords.append([atom1.x,atom1.y,atom1.z])
    secondcoords = []
    secondcoords.append([atom2.x,atom2.y,atom2.z])
    firstcoords = np.array(firstcoords)
    secondcoords = np.array(secondcoords)

    dist = np.linalg.norm(firstcoords-secondcoords)
    return dist

class Residue:
    def __init__(self,num,chain,name,atoms):
        self.num = num
        self.chain = chain
        self.name = name
        self.atoms = atoms

    #function assumes only 1 calpha
    def ca(self):
        for atom in self.atoms:
            if atom.atomid == ' CA ':
                return atom
        print 'residue',self.name,self.num,'has no CA. exiting'
        exit()

class Atom:
    def __init__(self,record,num,atomid,ali,Acode,x,y,z,occupancy,tempfact,segid,element,charge):
        self.record = record
        self.num = num
        self.atomid = atomid
        self.ali = ali
        self.Acode = Acode
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.tempfact = tempfact
        self.segid = segid
        self.element = element
        self.charge = charge
