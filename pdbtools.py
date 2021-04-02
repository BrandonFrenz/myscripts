#!/usr/bin
import re
import amino_acids
import os
import numpy as np
import gzip
import string


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

def mutate_residue(residue, newname):
    backbones = get_backbones(residue)
    residue.atoms = backbones
    residue.name = newname
    return residue

def convert_to_gly(residue):
    mutate_residue(residue, 'GLY')
    return residue
    

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

def write_resis_to_pdb(resis,name,sort=True):
    pdblines = make_pdblines_from_residues(resis,sort)
    write_pdb(pdblines,name)

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
    y = float(''.join(ll[38:46]))
    z = float(''.join(ll[46:54]))
    occupancy = ''.join(ll[55:60])
    temp = ''.join(ll[60:66])
    segid = ''.join(ll[72:76])
    element = ''.join(ll[76:78])
    charge = ''.join(ll[79:80])
    atom = Atom(record,num,atomid,ali,achar,x,y,z,occupancy,temp,segid,element,charge)
    return atom

def get_unopened_residue_list(pdbfile):
    if pdbfile.endswith('.gz'):
        gfile = gzip.open(pdbfile)
        residues = get_residue_list(gfile.readlines())
        gfile.close()
    else:
        with open(pdbfile,'r') as pfile:
            residues = get_residue_list(pfile.readlines())
    return residues

def get_residue_list(pdbfile):
    residues = []
    previousresid = ('','')
    previousicode = ''
    previousname = ''
    currentres = ''
    linecount = 0
    for line in pdbfile:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            ll = list(line)
            resnum = int(''.join(ll[22:26]))
            icode = ll[26]
            chainid = ''.join(ll[21])
            currentresid = (resnum,chainid)
            resname = ''.join(ll[17:20])
            if previousresid != currentresid or previousicode != icode or previousname != resname:   
                previousresid = currentresid
                previousicode = icode
                previousname = resname
                if currentres != '':
                    isterm = False
                    if pdbfile[linecount-1].startswith('TER'):
                        isterm = True
                    currentres.isterm = isterm
                    residues.append(currentres)
                atom = atom_from_pdbline(line)
                currentres = Residue(resnum,chainid,resname,icode,[])
                currentres.atoms.append(atom)
            else:
                atom = atom_from_pdbline(line)
                currentres.atoms.append(atom)
        linecount+=1
    if currentres != '':
        currentres.isterm = True
        residues.append(currentres)
    else:
        print('Warning: no residues found in', pdbfile)
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

#converts all modified residues to the standard format
def convert_to_standard_aas(residues):
    for res in residues:
        if res.name in amino_acids.modres.keys():
            res.name = amino_acids.modres[res.name]
            newatoms = []
            for atom in res.atoms:
                atom.record = 'ATOM  '
                newatoms.append(atom)
            res.atoms = newatoms

def strip_non_protein(residues):
    stripped = []
    for residue in residues:
        if residue.name in amino_acids.longer_names:
            stripped.append(residue)
    return stripped

def is_protein(residue):
    if residue.name in amino_acids.longer_names:
        return True
    else:
        return False


def set_full_occupancy(residues):
    for res in residues:
        for atom in res.atoms:
            atom.occupancy = " 1.00"

def clear_segid(residues):
    for res in residues:
        for atom in res.atoms:
            atom.segid = "    "

#this function will result in ter statements being added whereever non continuous residue numbering occurs
def add_ters_to_noncontres(residues):
    updated_residues = []
    for it in range(0, len(residues)-1):
        res = residues[it]
        if res.num != residues[it+1].num-1:
            res.isterm = True
        updated_residues.append(res)
    return updated_residues

def make_pdblines_from_residues(reslist,sort=True):
    pdbfile = []
    previousres = (0,'0')
    if sort:
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
    coord = float("{0:.3f}".format(coord))
    coord = str(coord).split('.')
    predec = coord[0]
    postdec = coord[1]
    newstr = ' '*(4-len(predec))+predec+'.'+postdec+'0'*(3-len(postdec))
    return newstr

def atomlist_rms(atoms1,atoms2):
    assert(len(atoms1) == len(atoms2)), 'unequal number of atoms exiting'
    n = len(atoms1)
    atomit = 0
    distsqr = 0
    while atomit < len(atoms1):
        dist = atom_dist(atoms1[atomit],atoms2[atomit])
        distsqr+=dist*dist
        atomit+=1
    rms = np.sqrt(distsqr/atomit)
    return rms

def atomlist_GDTha(atoms1,atoms2):
    assert(len(atoms1) == len(atoms2)), 'unequal number of atoms exiting'
    atomit = 0
    GDTha = 0.0
    while atomit < len(atoms1):
        dist = atom_dist(atoms1[atomit],atoms2[atomit])
        if dist < 0.5:
            GDTha+=1
        if dist < 1.0:
            GDTha+=1
        if dist < 2.0:
            GDTha +=1
        if dist < 4.0:
            GDTha +=1
        atomit+=1
    return GDTha/(atomit*4)

def get_backbones(res, include_cb=True):
    backbones = [' N  ',' CA ',' C  ', ' O  ']
    if include_cb:
        backbones.append(' CB ')
    resbb = []
    for atom in res.atoms:
        if atom.atomid in backbones:
            resbb.append(atom)
    return resbb

def atom_dist(atom1,atom2):
    firstcoords = []
    firstcoords.append([atom1.x,atom1.y,atom1.z])
    secondcoords = []
    secondcoords.append([atom2.x,atom2.y,atom2.z])
    firstcoords = np.array(firstcoords)
    secondcoords = np.array(secondcoords)

    dist = np.linalg.norm(firstcoords-secondcoords)
    return dist

def res_distance(res1, res2):
    closest_dist = None
    for atom1 in res1.atoms:
        for atom2 in res2.atoms:
            dist = atom_dist(atom1, atom2)
            if closest_dist == None or dist < closest_dist:
                closest_dist = dist
    return closest_dist

def connection_distance(res1, res2):
    res1C = None
    for atom in res1.atoms:
        if atom.atomid == ' C  ':
            res1C = atom
    res2N = None
    for atom in res2.atoms:
        if atom.atomid == ' N  ':
            res2N = atom
    return atom_dist(res1C, res2N)

def res_ca_distance(res1, res2):
    res1CA = None
    res2CA = None
    for atom in res1.atoms:
        if atom.atomid == ' CA ':
            res1CA = atom
    for atom in res2.atoms:
        if atom.atomid == ' CA ':
            res2CA = atom
    if res1CA == None or res2CA == None:
        print('No CA atoms found in one of the residues exiting', res1.num, res1.name, res2.num, res2.name)
        exit()
    else:
        return atom_dist(res1CA, res2CA)


def get_cas(residues):
    cas = []
    for res in residues:
        for atom in res.atoms:
            if atom.atomid == ' CA ':
                cas.append(atom)
    return cas

def get_ca(res):
    for atom in res.atoms:
        if atom.atomid == ' CA ':
            return atom

#returns a list containing the sequence of each chain
def get_sequences(residues):
    sequences = []
    current_sequence = []
    previouschain = ''
    for residue in residues:
        if residue.chain != previouschain and previouschain != '':
            sequences.append(current_sequence)
            current_sequence = []
        if residue.name in amino_acids.longer_names:
            current_sequence.append(amino_acids.longer_names[residue.name])
        else:
            current_sequence.append('X')
        previouschain = residue.chain
    if len(current_sequence) > 1:
        sequences.append(current_sequence)
    return sequences
            

def get_chain_resis(residues,chain):
    chain_resis = []
    for residue in residues:
        if residue.chain == chain:
            chain_resis.append(residue)
    return chain_resis

def get_chains(residues):
    chains = []
    for residue in residues:
        if residue.chain not in chains:
            chains.append(residue.chain)
    return chains

def relabel_chains(resis):
    relabeled = []
    chainindex = 0
    previous_chain = None
    for res in resis:
        if previous_chain != None and previous_chain != res.chain:
            chainindex+=1
        print(chainindex)
        newchain = string.ascii_uppercase[chainindex]
        previous_chain = res.chain
        res.chain = newchain
        relabeled.append(res)
    return relabeled

def write_fragments(fragments,filename):
    with open(filename,'w') as outfile:
        for fragment in fragments:
            if fragment.fragid == 'position:':
                outfile.write('position:'+' '*(13-len(str(fragment.num)))+str(fragment.num)+' '+'neighbors:'+' '*(13-len(str(fragment.neighbors)))+str(fragment.neighbors)+'\n')
                for line in fragment.lines:
                    outfile.write(line)

def parse_fragments(fragfile):
    if 'position' in fragfile[0]:
        fragments = parse_position_fragments(fragfile)
        return fragments

def parse_position_fragments(fragfile):
    full_fragments = []
    lines = []
    fragid = ''
    fragnum = ''
    neighbors = ''
    for line in fragfile:
        if line.startswith('position:'):
            if fragid != '':
                fragments = FRAGMENT_LIST(fragid,fragnum,neighbors,lines)
                full_fragments.append(fragments)
                lines = []
            line = line.split()
            fragid = 'position:'
            fragnum = int(line[1])
            neighbors = int(line[3])
        else:
            lines.append(line)

    fragments = FRAGMENT_LIST(fragid,fragnum,neighbors,lines)
    full_fragments.append(fragments)
    return full_fragments

def get_ca(resi):
    for atom in resi.atoms:
        if atom.atomid == ' CA ':
            return atom

def clean_pdbs(pdbs):
    for pdb in pdbs:
        residues = get_unopened_residue_list(pdb)
        clean_resis = []
        for residue in residues:
            if residue.name not in amino_acids.longer_names:
                continue
            else:
                clean_resis.append(residue)
        write_resis_to_pdb(clean_resis,pdb)

def get_sequence(residues, term_to_slash = True):
    sequence = []
    for residue in residues:
        if residue.name in amino_acids.longer_names:
            sequence.append(amino_acids.longer_names[residue.name])
        else:
            sequence.append('X')
        if residue.isterm and term_to_slash:
            sequence.append('/')
    return sequence

def make_point_to_MG(point,resnum):
    atom = Atom('ATOM  ',resnum,' MG ',' ',' ',point[0],point[1],point[2],' ',' ',' ','MG',' ')
    residue = Residue(resnum,'A','UNK',[atom])
    return residue

def get_resid(residue):
    return (residue.name,residue.chain,residue.num,residue.icode)

def backbone_rmsd(res1, res2, include_cb=True):
    bb1 = get_backbones(res1, include_cb)
    bb2 = get_backbones(res2, include_cb)
    if len(bb1) != len(bb2):
        print('warning backbones are not equal in length returning -1 for rmsd')
        return -1
    return atomlist_rms(bb1,bb2)

def get_residues_backbone_rmsd(reslist1, reslist2, include_cb=True):
    backbones1 = []
    backbones2 = []
    for res in reslist1:
        backbones1+=get_backbones(res, include_cb)
    for res in reslist2:
        backbones2+=get_backbones(res, include_cb)
    if len(backbones1) != len(backbones2):
        print('warning backbones are not equal in length returning -1 for rmsd')
        print(len(backbones1), len(backbones2))
        return -1
    return atomlist_rms(backbones1,backbones2)

#currently a very crude representation of all the fragments
class FRAGMENT_LIST:
    #fragnum corresponds to the first residue in the fragment
    def __init__(self,fragid,num,neighbors,lines):
        self.fragid = fragid
        self.num = num
        self.neighbors = neighbors
        self.lines = lines

class Residue:
    def __init__(self,num,chain,name,icode,atoms):
        self.num = num
        self.chain = chain
        self.name = name
        self.icode = icode
        self.atoms = atoms

    #function assumes only 1 calpha
    def ca(self):
        for atom in self.atoms:
            if atom.atomid == ' CA ':
                return atom
        print('residue',self.name,self.num,'has no CA. exiting')
        exit()

    def is_het(self):
        if self.atoms[0].record == 'HETATM':
            return True
        else:
            return False

    #def __eq__(self,other):
    #    if type(other) == type(self):
    #        return self.__dict__ == other.__dict__#self.num == other.num and self.name == other.name and self.chain == other.chain and self.icode == other.icode
    #    return False

    #def __ne__(self,other):
    #    return not self.__eq__(other)

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
