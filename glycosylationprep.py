#!/usr/bin/python
import argparse
import pdbtools
import copy
import amino_acids
import operator
import re
import numpy
from collections import defaultdict

def main():
    args = parseargs()
    if args.residues is not None:
        args.residues = parse_residues(args)
    if args.method == 'asym':
        if args.chains == 'ALL':
            print "you need to select the chains of the asymmetric unit with -c"
            exit()
        grab_asym_unit(args)
    if args.method == 'renum':
        renumber_amino_acids(args)
    if args.method == 'subsugs':
        substitute_sugars(args)
    if args.method == 'copy sugars':
        copy_sugars(args)
    if args.method == 'stripsug':
        strip_non_aa(args)
    if args.method == 'striph':
        strip_hydrogens(args)
    if args.method == 'collect':
        collect_sugars(args)
    if args.method == 'fixlinks':
        fix_all_links(args)
    if args.method == 'ters':
        add_ters(args)
    if args.method == 'individualres':
        renum_individual_res(args)
    if args.method == 'verify':
        verify_links(args)
    if args.method == 'convert':
        convert_types(args)
    if args.method == 'indisub':
        substitute_individual_resis(args)
    if args.method == 'remove':
        remove_residues(args)
    if args.method == 'add':
        add_indi_res(args)
    if args.method == 'rnum':
        print_rosnum(args)
    if args.method == 'pnum':
        print_pdbnum(args)
    if args.method == 'filterlinks':
        filter_links(args)
    if args.method == 'stripanomer':
        remove_anomeric_atoms(args)
    if args.method == 'atm':
        convert_atomrecord(args)
    if args.method == 'remove_dup_links':
        remove_duplicate_links(args)
    if args.method == 'checklinks':
        check_links(args)
    if args.method == 'fix chain numbers':
        fix_chain_numbers(args)
    if args.method == 'reorder residues':
        reorder_residues(args)
    if args.method == 'grab residues':
        grab_residues(args)
    if args.method == 'match nums':
        match_nums(args)
        

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb',help='The pdb to resolve')
    parser.add_argument('-pp','--positionpdb',help='The pdb to steal the residdues from')
    parser.add_argument('-o','--output',default='prepped_glysoylations.pdb',help='The name of your output')
    parser.add_argument('-a','--add',type=int,default=0,help='The number to add/subtract to every residue')
    parser.add_argument('-c','--chains',nargs="+",default='ALL',help='The chains to use')
    parser.add_argument('-m','--method',help='What method would you like to use. Read script "main" function for full list')
    parser.add_argument('-sp','--sugarpdbs',nargs="+",help='Collects all the non AA residues and their link records and attaches them to the first pdb')
    parser.add_argument('-r','--residues',nargs="+",help='The residues to change')
    parser.add_argument('-nr','--newresnums',nargs="+",type=int,help='The matching number for the new residues')
    parser.add_argument('-f','--formats',help='The type of format you want, rosetta or pdb')
    parser.add_argument('-dc','--distancecut',type=float,default=3.5,help='The distance cutoff for links')
    parser.add_argument('-orl','--override_links',type=int,default=0,help='Over ride the error that prevents deletion of links across chains being removed')
    args = parser.parse_args()
    return args

def parse_residues(args):
    residues = []
    for res in args.residues:
        res = res.split('-')
        chain = re.split('\d+',res[0])[0]
        lower = int(re.split('(\d+)',res[0])[1])
        upper = lower
        if len(res) > 1:
            upper = int(res[1])
        for i in range(lower,upper+1):
                residues.append(i)
    return residues

#catch all for anything in the header no included other places
def get_header(pdbfile):
    header = True
    headerlines = []
    for line in pdbfile:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            header = False
        linkline = False
        if line.startswith('LINK') or line.startswith('HETNAM') or line.startswith('SSBOND'):
            linkline = True
        if header == True and not linkline:
            headerlines.append(line)
    return headerlines


def get_tailer(pdbfile):
    intail = False
    lastatom = 0
    it = 0
    for line in pdbfile:
        if line.startswith('ATOM') or line.startswith('HETATM'):
                lastatom = it
        it+=1
    tail = []
    tailcount = lastatom+1
    while tailcount < len(pdbfile):
        if pdbfile[tailcount].startswith('TER'):
            tailcount+=1
            continue
        tail.append(pdbfile[tailcount])
        tailcount+=1
    return tail

def read_links(pdbfile):
    links = []
    for line in pdbfile:
        if line.startswith('LINK') or line.startswith('HETNAM'):
            try:
                link = create_link_from_pdbline(line)
                if 'LINK' in link.record:
                    links.append(link)
                elif link.record == 'HETNAM':
                    links.append(link)
            except:
                print 'skipping possible link',line
                continue
    return links

#creates a hetnam or link object based on line
def create_link_from_pdbline(line):
    line = line.rstrip('\n')
    if line.startswith('LINK'):
        ll = list(line)
        record = ''.join(ll[0:6])
        latom = ''.join(ll[12:16])
        lali = ll[16]
        lres = ''.join(ll[17:20])
        lchain = ll[21]
        lnum = int(''.join(ll[22:26]))
        licode = ll[26]
        uatom = ''.join(ll[42:46])
        uali = ll[46]
        ures = ''.join(ll[47:50])
        uchain = ll[51]
        unum = int(''.join(ll[52:56]))
        uicode = ll[56]
        sym1 = ''.join(ll[59:65])
        sym2 = ''.join(ll[66:72])
        dist = ''.join(ll[74:78])
        link = LINK(record,latom,lali,lres,lchain,lnum,licode,uatom,uali,ures,uchain,unum,uicode,sym1,sym2,dist)
        return link
    elif line.startswith('HETNAM'):
        line = line.strip()
        ll = list(line)
        record = ''.join(ll[0:6])
        res = ''.join(ll[11:14])
        chain = ''.join(ll[15])
        num = int(''.join(ll[16:20]))
        linkage = ''.join(ll[22:39])
        link = HETNAM(record,res,chain,num,linkage)
        return link
    else:
        print "You have attempted to create a link object from an unsupported pdbline " + str(line) + " exiting script"
        exit()

def resolve_backward_links(links):
    fixedlinks = []
    for link in links:
        if link.lnum > link.unum:
            link = invert_link(link)
        fixedlinks.append(link)
    return fixedlinks

def invert_link(link):
    newlink = copy.deepcopy(link)
    newlink.latom = link.uatom
    newlink.lali = link.uali
    newlink.lres = link.ures
    newlink.lchain = link.uchain
    newlink.lnum = link.unum
    newlink.licode = link.uicode
    newlink.uatom = link.latom
    newlink.uali = link.lali
    newlink.ures = link.lres
    newlink.uchain = link.lchain
    newlink.unum = link.lnum
    newlink.uicode = link.licode
    return newlink

def renumber_residue_and_corresponding_links(residue,newnum,links,disulfides):
    newres = copy.deepcopy(residue)
    newres.num = newnum
    newlinks = []
    for link in links:
        if 'LINK' in link.record:
            if link.lnum == residue.num:
                link.lnum = newres.num
            if link.unum == residue.num:
                link.unum = newres.num
        elif 'HETNAM' in link.record:
            if link.num == residue.num:
                link.num = newres.num
        newlinks.append(link)
    new_disulfides = []
    for disulf in disulfides:
        if disulf.lnum == residue.num:
            disulf.lnum = newres.num
        if disulf.unum == residue.num:
            disulf.unum = newres.num
        new_disulfides.append(disulf)
    return newres,newlinks,new_disulfides

#this function will renumber all the amino acids residues and their corresponding link and disulfide records
def renumber_amino_acids(args):
    #read pdb
    header,links,disulfides,residues,tailer = parse_pdbfile(args.pdb)

    #renumber residues
    newres = []
    for res in residues:
        if res.name in amino_acids.longer_names:
            res,links,disulfides = renumber_residue_and_corresponding_links(res,res.num+args.add,links,disulfides)
            newres.append(res)
        else:
            newres.append(res)
    residues = newres

    write_pdb(args,header,links,disulfides,residues,tailer)

#takes a LINK object and returns a string suitable for a pdbfile
def create_pdbline_from_link(link):
    if 'LINK' in link.record:
        line = list(' '*80)
        line[0:6] = link.record
        line[12:16] = link.latom
        line[16] = link.lali
        line[17:20] = link.lres
        line[21] = link.lchain
        line[22:26] = ' '*(4-len(str(link.lnum)))+str(link.lnum)
        line[26] = link.licode
        line[42:46] = link.uatom
        line[46] = link.uali
        line[47:50] = link.ures
        line[51] = link.uchain
        line[52:56] = ' '*(4-len(str(link.unum)))+str(link.unum)
        line[56] = link.uicode
        line[59:65] = link.sym1
        line[66:72] = link.sym2
        line[73:78] = link.dist.rstrip('\n')
        line+='\n'
        return ''.join(line)
    elif link.record == 'HETNAM':
        line = list(' '*80)+list('\n')
        line[0:6] = link.record
        line[11:14] = link.res
        line[15] = link.lchain
        line[16:20] = ' '*(4-len(str(link.num)))+str(link.num)
        line[22:38] = link.linkage
        return ''.join(line)
    else:
        print link.record

def fix_all_links(args):

    header,links,disulfs,residues,tailer  = parse_pdbfile(args.pdb)
    
    links = resolve_backward_links(links)
    links = reorder_links(links)
    write_pdb(args,header,links,disulfs,residues,tailer)

def reorder_links(links):
    newlinks = []
    while len(newlinks) < len(links):
        current_lowest = None
        for link in links:
            if link in newlinks:
                continue
            if current_lowest == None:
                current_lowest = link
            if current_lowest.lchain > link.lchain:
                current_lowest = link
            elif current_lowest.lnum > link.lnum:
                current_lowest = link
            elif current_lowest.lnum == link.lnum:
                if current_lowest.unum > link.unum:
                    current_lowest = link
        print current_lowest.lchain
        newlinks.append(current_lowest)
    return newlinks
    return sorted(links, key = lambda x: (x.lchain, x.lnum, x.unum))

def substitute_sugars(args):
    #positions
    positionpdb = open(args.positionpdb,'r').readlines()
    position_residues = pdbtools.get_residue_list(positionpdb)
    #sugars
    sugar_header,sugar_links,sugar_disulfides,sugar_residues,sugar_tail = parse_pdbfile(args.pdb)
    #add residues not accounted for
    for res in sugar_residues:
        hasalready = False
        #don't add protein residues
        if res.name in amino_acids.longer_names:
            hasalready = True
        for posres in position_residues:
            if posres.num == res.num and posres.chain == res.chain and posres.icode == res.icode:
                hasalready = True
        if hasalready == False:
            position_residues.append(res)
    if 'ALL' not in args.chains:
        only_chain_res = []
        for res in position_residues:
            if res.chain in args.chains:
                only_chain_res.append(res)
        position_residues = only_chain_res

    write_pdb(args,sugar_header,sugar_links,sugar_disulfides,position_residues,sugar_tail)

#rather than using the coordinates this function simple pastes the glycans blind into the new pdb
def copy_sugars(args):
    positionpdb = pdbtools.get_unopened_residue_list(args.positionpdb)
    sugar_header,sugar_links,sugar_disulfides,sugar_residues,sugar_tail = parse_pdbfile(args.pdb)
    nonaas = []
    for res in sugar_residues:
        if res.name not in amino_acids.longer_names:
            nonaas.append(res)
    newpdb = positionpdb+nonaas
    pdbtools.write_resis_to_pdb(newpdb,args.output)

#renumber the residues to match the positionpdb numbers
def match_nums(args):
    positionpdb = pdbtools.get_unopened_residue_list(args.positionpdb)
    sugar_header,sugar_links,sugar_disulfides,sugar_residues,sugar_tail = parse_pdbfile(args.pdb)
    nonaas = []
    it = 0
    print len(positionpdb),len(sugar_residues)
    while it <len(sugar_residues):
        #sugar_residues[it].num = positionpdb[it].num
        sugar_residues[it],sugar_links,sugar_disulfides = renumber_residue_and_corresponding_links(sugar_residues[it],positionpdb[it].num,sugar_links,sugar_disulfides)
        it+=1
    write_pdb(args,sugar_header,sugar_links,sugar_disulfides,sugar_residues,sugar_tail)

def parse_pdbfile(pdb):
    pdbfile = open(pdb,'r').readlines()
    links = read_links(pdbfile)
    header = get_header(pdbfile)
    tailer = get_tailer(pdbfile)
    residues = pdbtools.get_residue_list(pdbfile)
    disulfs = get_disulfides_from_pdb(pdbfile)
    return header,links,disulfs,residues,tailer

def get_disulfides_from_pdb(pdbfile):
    disulfides = []
    for line in pdbfile:
        if line.startswith('SSBOND'):
            line = line.strip()
            disulf = get_disulfide_from_line(line)
            disulfides.append(disulf)
    return disulfides

def get_disulfide_from_line(line):
        ll = list(line)
        record = ''.join(ll[0:6])
        bondnum = ''.join(ll[7:10])
        lchain = ll[15]
        lresnum = int(''.join(ll[17:21]).strip())
        licode = ll[21]
        uchain = ll[29]
        uresnum = int(''.join(ll[31:35]).strip())
        uicode = ll[35]
        symm1 = ll[59:65]
        symm2 = ll[66:72]
        try:
            length = float(''.join(ll[74:78]))
        except:
            length = ''.join(ll[74:78])
        disulf = DISULFIDE(record,bondnum,lchain,lresnum,licode,uchain,uresnum,uicode,symm1,symm2,length)
        return disulf

def get_lines_from_disulfides(disulfides):
    lines = []
    for disulfide in disulfides:
        line = make_line_from_disulf(disulfide)
        lines.append(line)
    return lines

def make_line_from_disulf(disulf):
    line = list(' '*78)
    line[0:6] = disulf.record
    line[7:10] = disulf.bondnum
    line[11:14] = ['C','Y','S']
    line[15] = disulf.lchain
    line[17:21] = ' '*(4-len(str(disulf.lnum)))+str(disulf.lnum)
    line[21] = disulf.licode
    line[25:28] = ['C','Y','S']
    line[29] = disulf.uchain
    line[31:35] = ' '*(4-len(str(disulf.unum)))+str(disulf.unum)
    line[35] = disulf.uicode
    line[59:65] = disulf.symm1
    line[66:72] = disulf.symm2
    line[73:78] = ' '*(4-len(str(disulf.length)))+str(disulf.length)
    line = ''.join(line)+"\n"
    return line

def grab_asym_unit(args):
    #read parse the pdbfile
    header,links,disulfides,residues,tailer = parse_pdbfile(args.pdb)
    
    #filter asymm residues
    newresidues = []
    for res in residues:
        if res.chain not in args.chains:
            continue
        newresidues.append(res)
    residues = newresidues

    #filter asym links
    asymlinks = []
    for link in links:
        if link.record == 'HETNAM':
            if link.lchain in args.chains:
                asymlinks.append(link)
        if 'LINK' in link.record:
            if link.lchain in args.chains and link.uchain in args.chains:
                asymlinks.append(link)
            elif link.lchain in args.chains or link.uchain in args.chains:
                print "Invalid attempt to create a link between a residue on a chain you intended to keep and one being removed"
                print link.lchain,link.lnum,link.uchain,link.unum
                if args.override_links != 1:
                    exit()
            else:
                continue
    links = asymlinks
 
 #filter asym disulfides
    asymdisulf = []
    for disulf in disulfides:
        if disulf.lchain and disulf.uchain in args.chains:
            asymdisulf.append(disulf)
        elif disulf.lchain in args.chains or disulf.uchain in args.chains:
            print "Invalid attemp to create disulfide bond between a residue on a chain you intend to keep and one being removed"
            print disulf.lchain,disulf.lnum,disulf.uchain.disulf.unum
    disulfides = asymdisulf


    write_pdb(args,header,links,disulfides,residues,tailer)

def write_pdb(args,header,links,disulfides,residues,tailer):
    #get the lines
    linklines = get_lines_from_links(links)
    reslines = pdbtools.make_pdblines_from_residues(residues)
    disulfides = get_lines_from_disulfides(disulfides)
    #write pdb
    newpdbfile = header+linklines+disulfides+reslines+tailer
    with open(args.output,'w') as outfile:
        for line in newpdbfile:
            outfile.write(line)

def get_lines_from_links(links):
    lines = []
    #links = sorted(links, key=operator.itemgetter(links.record,links.lnum))
    #links = sorted(links, key=lambda x: (x.record,x.lnum))
    for link in links:
        lines.append(create_pdbline_from_link(link))
    return lines

def strip_non_aa(args):
    header,links,disulfides,residues,tailer = parse_pdbfile(args.pdb)
    newres = []
    for res in residues:
        if res.name not in amino_acids.longer_names:
            continue
        else:
            newres.append(res)
    residues = newres
    pdbfile = pdbtools.make_pdblines_from_residues(residues)
    pdbtools.write_pdb(pdbfile,args.output)

def collect_sugars(args):
    sugars = []
    suglinks = []
    for pdb in args.sugarpdbs:
        header,links,disulfs,residues,tailer = parse_pdbfile(pdb)
        for res in residues:
            if res.name not in amino_acids.longer_names:
                has_sug = False
                for sug in sugars:
                    if sug.num == res.num:
                        has_sug = True
                if not has_sug:
                    sugars.append(res)
                    for link in links:
                        if link.lnum == res.num or link.unum == res.num:
                            suglinks.append(link)
    pdb1 = open(args.sugarpdbs[0]).readlines()
    residues = pdbtools.get_residue_list(pdb1)
    canonical_aas = []
    for res in residues:
        if res.name in amino_acids.longer_names:
            canonical_aas.append(res)
    residues = canonical_aas
    residues+=sugars
    suglinks = resolve_backward_links(suglinks)
    write_pdb(args,[],suglinks,[],residues,[])
        
def add_ters(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    it = 0
    ter_res = []
    while it< len(residues)-2:
        if residues[it].num != residues[it+1].num-1 and residues[it].num != residues[it+1].num:
            residues[it].isterm = True
        it+=1
    write_pdb(args,header,links,disulfs,residues,tailer)

def renum_individual_res(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    renumed = []
    for res in residues:
        if res.num in args.residues and (res.chain in args.chains or args.chains == 'ALL'):
            resindex = args.residues.index(res.num)
            newres = args.newresnums[resindex]
            res,links,disulfs = renumber_residue_and_corresponding_links(res,newres,links,disulfs)
        renumed.append(res)
    residues = renumed
    write_pdb(args,header,links,disulfs,residues,tailer)

def verify_links(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    for link in links:
        has_lower = False
        has_upper = False
        for res in residues:
            if res.num == link.lnum:
                has_lower = True
            if res.num == link.unum:
                has_upper = True
        if has_lower == False:
            print 'link',link.lnum,'-',link.unum,'has no lower corresponding residue'
        if has_upper == False:
            print 'link',link.lnum,'-',link.unum,'has no upper corresponding residue'
        if has_upper and has_lower:
            print link.lnum,'is valid'
    for link in links:
        linkcounter = 0
        for slink in links:
            if link.lnum == slink.lnum:
                linkcounter+=1
        if linkcounter > 1:
            print link.lnum,'is beta branched'

def convert_types(args):
    print "this is a one off function meant for a specific use case, it may not be applicable to you"
    if args.formats == None:
        print 'no format selected exiting protocol, use "-f rosetta/pdb"'
        exit()
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    newresidues = []
    resiudes = convert_MAN_to_pdb_format(residues,links)
    resit=-1
    for res in residues:
        resit+=1
        if args.formats == 'pdb':
            res.name = resid_to_pdbformat(res.name)
            if len(residues) > resit+1:
                if res.chain != residues[resit+1].chain:
                    res.isterm = True
                else:
                    res.isterm = False
            else:
                res.isterm = True
        newatoms = []
        for atom in res.atoms:
            if args.formats == 'rosetta':
                if atom.atomid.strip() == 'C7':
                    atom.atomid = ' CN2'
                if atom.atomid.strip() == 'C8':
                    atom.atomid = 'CAN2'
                if atom.atomid.strip() == 'O7':
                    atom.atomid = 'OCN2'
            if args.formats == 'pdb':
                atom.record == 'HETATM'
                if atom.atomid.strip() == 'CN2':
                    atom.atomid = ' C7 '
                if atom.atomid.strip() == 'CAN2':
                    atom.atomid = ' C8 '
                if atom.atomid.strip() == 'OCN2':
                    atom.atomid = ' O7 '
            if atom.element == ' H':
                continue
            newatoms.append(atom)
        res.atoms = newatoms
        newresidues.append(res)
    newlinks = []
    for link in links:
        if args.formats == 'pdb':
            if link.record == 'HETNAM':
                newlink = convert_HETNAM_to_LINK(args,link)
                newlinks.append(newlink)
            elif link.record.strip() == 'LINK':
                link.lres = resid_to_pdbformat(link.lres)
                link.ures = resid_to_pdbformat(link.ures)
                newlinks.append(link)
    links = newlinks
    residues = newresidues
    write_pdb(args,header,links,disulfs,residues,tailer)

def convert_MAN_to_pdb_format(residues,links):
    newresidues = []
    for res in residues:
        if res.name == 'Man':
            for link in links:
                if link.lnum == res.num:
                    if 'alpha' in link.linkage:
                        res.name = 'MAN'
                        break
                    if 'beta' in link.linkage:
                        res.name = 'BMA'
                        break
        newresidues.append(res)
    return newresidues

def resid_to_pdbformat(name):
    if name == 'Man':
        return 'BMA'
    if name == 'Glc':
        return 'NAG'
    else:
        return name

def convert_HETNAM_to_LINK(args,link):
    #def __init__(self,record,res,chain,num,linkage):
    oxynum = re.split('>|\)',link.linkage)[1]
    oxy = ' O'+oxynum+''
    ltype = link.res
    utype = ''
    if args.method == 'convert' and args.formats == 'pdb':
        if 'Man' in link.linkage:
            utype = 'BMA'
        if 'Glc' in link.linkage:
            utype = 'NAG'
        if link.res == 'Glc':
            ltype = 'NAG'
        if link.res == 'Man':
            ltype = 'BMA'
    else:
        print 'this function currently only does switching to pdb residue types'
        exit()
    #def __init__(self,record,latom,lali,lres,lchain,lnum,licode,uatom,uali,ures,uchain,unum,uicode,sym1,sym2,dist):
    newlink = LINK('LINK  ',oxy,' ',ltype,link.chain,link.num,' ',' C1',' ',utype,link.chain,link.num+1,' ','','','')
    return newlink

def strip_hydrogens(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    newresidues = []
    for res in residues:
        newatoms = []
        for atom in res.atoms:
            if atom.element.strip() == 'H':
                continue
            newatoms.append(atom)
        res.atoms = newatoms
        newresidues.append(res)
    residues = newresidues
    write_pdb(args,header,links,disulfs,residues,tailer)

def substitute_individual_resis(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    subheader,sublinks,subdisulfs,subresidues,subtailer = parse_pdbfile(args.positionpdb)
    subres = []
    print args.residues
    for res in residues:
        if res.num in args.residues:
            for sres in subresidues:
                #print 'subreso', res.num
                if sres.num == res.num:
                    res = sres
        subres.append(res)
    residues = subres
    write_pdb(args,header,links,disulfs,residues,tailer)

def add_indi_res(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    subheader,sublinks,subdisulfs,subresidues,subtailer = parse_pdbfile(args.positionpdb)
    for res in subresidues:
        if res.num in args.residues:
            residues.append(res)
    write_pdb(args,header,links,disulfs,residues,tailer)

def remove_residues(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    keptres = []
    for res in residues:
        if res.num in args.residues:
            continue
        keptres.append(res)
    residues = keptres
    write_pdb(args,header,links,disulfs,residues,tailer)

def print_rosnum(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    rosettanum = 1
    for res in residues:
        if res.num in args.residues:
            print res.num,'is',rosettanum,'in rosetta numbering'
        rosettanum+=1

def print_pdbnum(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    rosettanum = 1
    for res in residues:
        if rosettanum in args.residues:
            print res.name,res.num,res.chain,'is',rosettanum,'in rosetta numbering'
        rosettanum+=1

def filter_links(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    newlinks = []
    for link in links:
        if link.record.strip() != 'LINK':
            continue
        latom = ''
        uatom = ''
        for res in residues:
            if link.lnum == res.num:
                for atom in res.atoms:
                    if atom.atomid == link.latom:
                        latom = atom
            if link.unum == res.num:
                for atom in res.atoms:
                    if atom.atomid == link.uatom:
                        uatom = atom
        if latom == '' or uatom == '':
            print 'skipping link',link.lnum
            continue
        dist = atom_dist(latom,uatom)
        if dist > args.distancecut:
            print 'link',link.lnum,'to',link.unum,'is too far apart removing it. Overwrite this behavior by setting the distance cutoff, -dc, arbitrarily large'
            continue
        newlinks.append(link)
    links = newlinks
    write_pdb(args,header,links,disulfs,residues,tailer)

def atom_dist(atom1,atom2):
    a1 = numpy.array((atom1.x,atom1.y,atom1.z))
    a2 = numpy.array((atom2.x,atom2.y,atom2.z))
    return numpy.linalg.norm(a1-a2)

def remove_anomeric_atoms(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    newres = []
    #badatoms = ['O2','O3','O4','O6','C6']
    badatoms = ['O6','C6']
    for res in residues:
        if res.num in args.residues:
            newatoms = []
            for atom in res.atoms:
                if atom.atomid.strip() in badatoms or atom.element == ' H':
                    continue
                newatoms.append(atom)
            res.atoms = newatoms
        newres.append(res)
    residues = newres
    write_pdb(args,header,links,disulfs,residues,tailer)

def convert_atomrecord(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    newres = []
    for res in residues:
        atoms = []
        for atom in res.atoms:
            atom.record = 'ATOM  '
            atoms.append(atom)
        res.atoms = atoms
        newres.append(res)
    residues = newres
    write_pdb(args,header,links,disulfs,residues,tailer)

def remove_duplicate_links(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    newlinks = []
    for link in links:
        nondup = True
        for nlink in newlinks:
            if link.lnum == nlink.lnum and link.lchain == nlink.lchain and link.unum == nlink.unum and link.uchain == nlink.uchain:
                nondup = False
                print link.lnum,'-',link.unum,'is duplicate'
        if nondup:
            newlinks.append(link)
    links = newlinks
    write_pdb(args,header,links,disulfs,residues,tailer)

def check_links(args):
    header,links,disulfs,residues,tailer = parse_pdbfile(args.pdb)
    for res in residues:
        if res.name not in amino_acids.longer_names:
            haslink = False
            print 'checking',res.num,res.name
            for link in links:
                if res.num == link.lnum or link.unum:
                    haslink = True
            if haslink == False:
                print res.num,'is missing links'

def parse_connect_records(pdbfile):
    atom_connections = defaultdict(list)
    for line in pdbfile:
        if line.startswith('CONECT'):
            line = line.strip()
            atomid = int(line[6:11])
            connect1 = line[11:16]
            connect2 = line[16:21]
            connect3 = line[21:26]
            connect4 = line[26:31]
            connections = [connect1,connect2,connect3,connect4]
            for connection in connections:
                try:
                    if int(connection) > atomid:
                        atom_connections[atomid].append(int(connection))
                except:
                    continue
    return atom_connections

def fix_chain_numbers(args):
    pdbfile = open(args.pdb).readlines()
    residues = pdbtools.get_unopened_residue_list(args.pdb)
    connections = parse_connect_records(pdbfile)
    res_connections = map_atom_connections_to_residue(residues,connections)
    reslocation = -3
    mid_chain = []
    for cres in (sorted(res_connections.keys(), key=operator.attrgetter('num'))):
        cresid = (cres.name,cres.num,cres.chain,cres.icode)
        if cresid in mid_chain:
            continue
        chains = trace_chain(cres,res_connections)
        for chain in chains:
            print 'new chain'
            for res in chain:
                resid = (res.name,res.num,res.chain,res.icode)
                print resid
                mid_chain.append(resid)

def trace_chain(res,res_connections):
    connected = [res]
    visited = {}
    path = []
    chains = []
    while len(connected) > 0:
        current_res = connected.pop()
        if len(res_connections[current_res]) == 0:
            path.append(current_res)
            chains.append(copy.deepcopy(path))
            path.pop()
            continue
        if current_res not in visited.keys():
            visited[current_res] = 1
        else:
            visited[current_res] += 1
        if len(res_connections[current_res]) == visited[current_res]-1:
            continue
        else:
            connected_res = res_connections[current_res][visited[current_res]-2]
            connected.append(connected_res)
        connected.append(current_res)
        if current_res not in path:
            path.append(current_res)
    return chains

def map_atom_connections_to_residue(residues,connections):
    connected_atoms = {}
    res_connections = defaultdict(list)
    for atomnum in connections.keys():
        connected_resis = find_atomnum_residues([atomnum]+connections[atomnum],residues)
        residue = connected_resis[atomnum]
        #res_connections[connected_resis[connection].num] = []
        for res in connected_resis.keys():
            #if connected_resis[res] == residue:
            if res_matches(residue,connected_resis[res]):
                continue
            if connected_resis[res] not in res_connections[residue]:
                res_connections[residue].append(connected_resis[res])
    return res_connections
       
def res_matches(res1,res2):
    if res1.num == res2.num and res1.chain == res2.chain and res1.icode == res2.icode and res1.name == res2.name:
        return True
    else:
        return False

def find_atomnum_residues(nums,residues):
    matching_residues = {}
    for res in residues:
        for atom in res.atoms:
            if atom.num in nums:
                matching_residues[atom.num] = res
    return matching_residues

def fix_connections(branches):
    needsfixing = False
    branchnums = []
    for i in range(0,len(branches)):
        branch = branches[i]
        bnum = []
        for res in branch:
            if res.chain != branch[0].chain:
                print 'error',res.name,res.num,res.chain,res.icode,'does not share a chain with the first residue in this branch'
                exit()
            if res.icode != ' ':
                needsfixing = True

def reorder_residues(args):
    header,links,disulfides,residues,tailer = parse_pdbfile(args.pdb)
    current_chain = ''
    done = False
    while not done:
        done = True
        for i in range(1,len(residues)):
            if residues[i].chain != residues[i-1].chain:
                continue
            if residues[i].num < residues[i-1].num:
                residues[i], residues[i-1] = residues[i-1],residues[i]
                done = False
            if i+1 < len(residues):
                if residues[i].chain != residues[i+1].chain:
                    residues[i].isterm = True
                else:
                    residues[i].isterm = False
    write_pdb(args,header,links,disulfides,residues,tailer)

#has no chain specifier super hacky
def grab_residues(args):
    header,links,disulfides,residues,tailer = parse_pdbfile(args.pdb)
    newresidues = []
    for res in residues:
        if res.num in args.residues and (res.chain in args.chains or args.chains == 'ALL'):
            newresidues.append(res)
    newlinks = []
    for link in links:
        if (link.lnum in args.residues or link.unum in args.residues) and (link.lchain in args.chains or args.chains == 'ALL'):
            print 'appending link'
            newlinks.append(link)
    newdisulfs = []
    for disulf in disulfides:
        if disulf.lnum in args.residues or disulf.unum in args.residues and (disulf.lchain in args.chains or args.chains  == 'ALL'):
            newdisulfs.append(disulf)
    write_pdb(args,header,newlinks,newdisulfs,newresidues,tailer)



class LINK:
    def __init__(self,record,latom,lali,lres,lchain,lnum,licode,uatom,uali,ures,uchain,unum,uicode,sym1,sym2,dist):
        self.record = record
        self.latom = latom
        self.lali = lali
        self.lres = lres
        self.lchain = lchain
        self.lnum = lnum
        self.licode = licode
        self.uatom = uatom
        self.uali = uali
        self.ures = ures
        self.uchain = uchain
        self.unum = unum
        self.uicode = uicode
        self.sym1 = sym1
        self.sym2 = sym2
        self.dist = dist

class DISULFIDE:
    def __init__(self,record,bondnum,lchain,lnum,licode,uchain,unum,uicode,symm1,symm2,length):
        self.record = record
        self.bondnum = bondnum
        self.lchain = lchain
        self.lnum = lnum
        self.licode = licode
        self.uchain = uchain
        self.unum = unum
        self.uicode = uicode
        self.symm1 = symm1
        self.symm2 = symm2
        self.length = length

class HETNAM:
    def __init__(self,record,res,chain,num,linkage):
        self.record = record
        self.res = res
        self.lchain = chain
        self.num = num
        self.linkage = linkage
        self.lnum = num
        self.unum = num+1

main()
