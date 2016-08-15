#!/usr/bin
import re
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB import Residue
from Bio.PDB import StructureAlignment
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import amino_acids
import os

def get_pdb_id(pdbpath):
    pdbid = re.split('/|.pdb',pdbpath)[-2]
    return pdbid

#Makes a file named alignment.fasta suitable for clustalw input
def make_alignment_input(seq1,seq2):
    newfile = []
    newfile.append("> seq1 " + str(len(seq1)) + " residues\n")
    newfile.append(seq1+"\n")
    newfile.append("> seq2 " + str(len(seq2)) + " residues\n")
    newfile.append(seq2+"\n")
    with open ('alignment.fasta','w') as af:
        af.write(''.join(newfile))

#Makes and biopython alignment from two pdbfile inputs
def make_structure_alignment(pdb1,pdb2):
    seq1 = str(get_sequence(pdb1))
    seq2 = str(get_sequence(pdb2))
    s1 = get_bio_structure(pdb1)
    pdbid = get_pdb_id(pdb1)
    s2 = get_bio_structure(pdb2)
    make_alignment_input(seq1,seq2)#produces an alignment file of the 2 sequences named alignment.fasta
    cline = ClustalwCommandline('clustalw',infile='alignment.fasta')
    cline()
    falign = AlignIO.read("alignment.aln", "clustal")
    salign = StructureAlignment(falign,s1,s2)
    os.system('rm alignment.aln alignment.fasta alignment.dnd')
    return salign

def get_bio_structure(pdb):
    pdbid = get_pdb_id(pdb)
    parser = PDBParser()
    s = parser.get_structure(pdbid,pdb)
    return s

#takes a pdbfile and returns the sequence
def get_sequence(pdb):
    s = get_bio_structure(pdb)
    ppbuilder = PPBuilder()
    model_nr = 1
    poly_peps = ppbuilder.build_peptides(s,model_nr)
    for polypep in poly_peps:
        return polypep.get_sequence()

#takes a biopython residue and returns the backbone heavy atoms with a different residue type
def convert_bio_res(res,newtype):
    #resid = res.id
    #resname = newtype
    resid = res.id
    resname = newtype
    ressegid = res.segid
    backbone = ['N','C','CA','O']
    if res.get_resname() != 'Gly' or newtype == 'Gly':
        backbone.append(['CB'])
    newres = Residue.Residue(resid,resname,ressegid)
    for atom in res:
        if atom.get_name() in backbone:
            newres.add(atom)
    return newres 
        
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


class Residue:
    def __init__(self,num,chain,name,atoms):
        self.num = num
        self.chain = chain
        self.name = name
        self.atoms = atoms

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
