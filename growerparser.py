#!/usr/bin/python2.7

def parse_unopened_beamfile(beamfile):
    with open(beamfile,'r') as bf:
        beams = parse_beamfile(bf.readlines())
    return beams

def parse_beamfile(beamfile):
    data = []
    for line in beamfile:
        line = line.split()
        data.extend(line)
    it = 0
    solutionset = []
    while it < len(data):
        total_res = int(data[it])
        it+=1
        total_lower = int(data[it])
        it+=1
        total_upper = int(data[it])
        it+=1
        step = int(data[it])
        it+=1
        lower = data[it]
        it+=1
        torsions = []
        for i in range(0,total_res):
            #loop over the 3 torsional angles
            for ii in range(0,3):
                torsions.append(float(data[it]))
                it+=1
            nchi = int(data[it])
            it+=1
            torsions.append(nchi)
            for ii in range(0,nchi):
                torsions.append(float(data[it]))
                it+=1
        
        bbresidues = []
        for i in range(0,total_res):
            residue = []
            atomcount = int(data[it])
            it+=1
            for ii in range(0,atomcount):
                #atomxyz = []
                #loop over the 3 coordinates
                atomx = float(data[it])
                it+=1
                atomy = float(data[it])
                it+=1
                atomz = float(data[it])
                it+=1
                #for i in range(0,3):
                #    atomxyz.append(float(data[it]))
                #    it+=1
                beam_atom = BEAM_ATOM(ii,atomx,atomy,atomz)
                residue.append(beam_atom)
            bbresidues.append(residue)
        
        sheets = []
        nsheets = int(data[it])
        it+=1
        bonus_score = float(data[it])
        it+=1
        for i in range(0,nsheets):
            baseres = int(data[it])
            it+=1
            jumpid = int(data[it])
            it+=1
            rotation = []
            translation = []
            for ii in range(0,9):
                rotation.append(float(data[it]))
                it+=1
            for ii in range(0,3):
                translation.append(float(data[it]))
                it+=1
            nsheetres = int(data[it])
            it+=1
            sheet_torsions = []
            for ii in range(0,nsheetres):
                for i in range(0,3):
                    sheet_torsions.append(float(data[it]))
                    it+=1
                nchis = int(data[it])
                it+=1
                for ii in range(0,nchi):
                    sheet_torsions.append(float(data[it]))
                sheet = SHEET(baseres,jumpid,rotation,translation,nsheetres,sheet_torsions)
                sheets.append(sheet)
        score = float(data[it])
        it+=1
        rms = float(data[it])
        it+=1
        gdt = float(data[it])
        it+=1
        beam = BEAM(total_res,total_lower,total_upper,step,lower,torsions,bbresidues,sheets,bonus_score,score,rms,gdt)
        solutionset.append(beam)
    return solutionset
        
                
            
class BEAM:
    def __init__(self,total_res,total_lower,total_upper,step,lower,torsions,bbresidues,sheets,sheetbonus,score,rms,gdt):
        self.total_res = total_res
        self.total_lower = total_lower
        self.total_upper = total_upper
        self.step = step
        self.lower = lower
        self.torsions = torsions
        self.bbresidues = bbresidues
        self.sheets = sheets
        self.sheetbonus = sheetbonus
        self.score = score
        self.rms = rms
        self.gdt = gdt

class SHEET:
    def __init__(self,baseres,jumplabel,rotation,translation,total_residues,torsions):
        self.baseres = baseres
        self.jumplabel = jumplabel
        self.rotation = rotation
        self.translation = translation
        self.total_residues = total_residues
        self.torsions = torsions

class BEAM_ATOM:
    def __init__(self,atomid,x,y,z):
        self.atomid = atomid
        self.x = x
        self.y = y
        self.z = z
