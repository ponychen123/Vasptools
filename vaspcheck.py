#!/usr/bin/python3
#a tool to warn about common errors in vasp
#history: 2008-08-29 Peter Larsson original version
#         2012-12-11 Peter Larsson Rewritten and improved
#         2019-07-29 ponychen rewrite and little improved
#         further future may imporoved, now set as a template for my future code 
#author:ponychen
#20190729
#email:18709821294@outlook.com

import os
import sys
import subprocess
import logging
from optparse import OptionParser
import re, collections
import math
import random

#Class for VASP input/output files

def cross(a, b):
    #get the cross product of vectors a and b
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

def dot(a, b):
    #get the dot product of vector a and b
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

class Poscar():
    def __init__(self, filename="POSCAR", use_direct=True):
        #read data from POSCAR
        poscarfile = open(filename, 'r')
        poscardata = poscarfile.readlines()
        poscarfile.close()
        #look for POTCAR file and construct atom label list
        atoms = []
        potcar_lines = subprocess.getoutput("grep TITEL POTCAR").split("\n")
        if len(potcar_lines) == 0:
            print("Cannot find atomic information in POTCAR, Missing POTCAR?")
            exit()
        for line in potcar_lines:
            words = line.split()
            assert words[0] == 'TITEL', "Attention, your POTCAR format is wrong!!!"

            #Note, we need to split _ to deal with eg "Li_sv"
            atoms.append(words[3].split("_")[0])

        #save the POSCAR header tag for later use when saving the file
        self.description = poscardata[0].strip()

        self.direct = use_direct

        latticeconstant = float(poscardata[1].strip())

        #check if volume scaling (2nd line is negtive number) is used
        volumescaling = False
        if latticeconstant < 0.0:
            volumescaling = True

        #read cell parameters
        self.a = list(map(float, poscardata[2].split()))
        self.b = list(map(float, poscardata[3].split()))
        self.c = list(map(float, poscardata[4].split()))

        if volumescaling:
            #lattice constant is really a volume scaling factor
            #calculate unscaled volume from determinant

            unscaledvolume = self.a[0]*self.b[1]*self.c[2]-self.a[0]*self.b[2]*self.c[1]+\
                    self.a[1]*self.b[2]*self.c[0]-self.a[1]*self.b[0]*self.c[2]+\
                    self.a[2]*self.b[0]*self.c[1]-self.a[2]*self.b[1]*self.c[0]
            scalingfactor = (-latticeconstant/unscaledvolume)**(1.0/3.0)

            #apply scaling to lattice constants
            for i in range(3):
                self.a[i] *= scalingfactor
                self.b[i] *= scalingfactor
                self.c[i] *= scalingfactor
        else:
            #apply lattice constants
            for i in  range(3):
                self.a[i] *= latticeconstant
                self.b[i] *= latticeconstant
                self.c[i] *= latticeconstant
        
        #get the constrain information in POSCAR
        if poscardata[7].upper()[0] == "S":
            self.constrain = True
        else:
            self.constrain = False

        #read atomic positions
        #the "atoms" array is an list like this "("Fe", {0.95, 0.05, 0.001})"
        #read atom counts
        atomlabels = []
        atomcounts = list(map(int, poscardata[6].split()))
        n_atoms = 0
        for i in range(len(atomcounts)):
            n_atoms += atomcounts[i]
            for j in range(atomcounts[i]):
                atomlabels.append(atoms[i])

        offset = 7
        #check for selective dynamics
        if poscardata[offset].upper()[0] == "S":
            offset += 1

        #check for direct coordinates
        direct_in_file = True
        if poscardata[offset].upper()[0] == "C":
            direct_in_file = False

        #scan atom positions
        offset += 1

        self.atoms = []
        for i in range(offset, offset+n_atoms):
            if poscardata[i] != "" and not ''.join(poscardata[i].split()).isalpha():
                thisatom = [0,0]
                thisatom[0] = atomlabels[i-offset]

                rvector = []
                abc = poscardata[i].split()

                if len(abc) > 0:
                    if use_direct:
                        #want ouput direct positions
                        if direct_in_file:
                            #Just save coords
                            for j in range(3):
                                rvector.append(float(abc[j]))
                        else:
                            #transfer cartesian coordination to direct coordination
                            r,s,t = [0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]
                            k,l,m = [0,0,0,0],[0,0,0,0],[0,0,0,0]
                            r[3],s[3],t[3] = self.a
                            r[4],s[4],t[4] = self.b
                            r[5],s[5],t[5] = self.c
                            A1=s[5]*r[3]-s[3]*r[5]
                            B1=t[5]*r[3]-t[3]*r[5]
                            A2=s[4]*r[3]-s[3]*r[4]
                            B2=t[4]*r[3]-t[3]*r[4]
                            Q=B2*A1-B1*A2
                            k[3]=(A2*r[5]-A1*r[4])/Q
                            l[3]=A1*r[3]/Q
                            m[3]=-A2*r[3]/Q
                            k[2]=(B2*r[5]-B1*r[4])/(-Q)
                            l[2]=B1*r[3]/(-Q)
                            m[2]=-B2*r[3]/(-Q)
                            k[1]=(1-s[3]*k[2]-t[3]*k[3])/r[3]
                            l[1]=-(s[3]*l[2]+t[3]*l[3])/r[3]
                            m[1]=-(s[3]*m[2]+t[3]*m[3])/r[3]
                            rvector.append(k[1]*float(abc[0])+l[1]*float(abc[1])+m[1]*float(abc[2]))
                            rvector.append(k[2]*float(abc[0])+l[2]*float(abc[1])+m[2]*float(abc[2]))
                            rvector.append(k[3]*float(abc[0])+l[3]*float(abc[1])+m[3]*float(abc[2]))
                    else:
                        if direct_in_file:
                            #convert direct coordinates to cartesian here
                            rvector.append(float(abc[0])*self.a[0]+float(abc[1])*self.a[1]+float(abc[2])*self.a[2])
                            rvector.append(float(abc[0])*self.b[0]+float(abc[1])*self.b[1]+float(abc[2])*self.b[2])
                            rvector.append(float(abc[0])*self.c[0]+float(abc[1])*self.c[1]+float(abc[2])*self.c[2])

                        else:
                            rvector.append(float(abc[0]))
                            rvector.append(float(abc[1]))
                            rvector.append(float(abc[2]))

                    thisatom[1] = rvector
                    self.atoms.append(thisatom)
            else:
                print("Reached end of file while read POSCAR " + str(i) + "th lines.\n")
        
        #get the fix information of the atoms
        if self.constrain:
            self.fix = []
            for i in range(9, 9+n_atoms):
                string = poscardata[i].split()[3:6]
                self.fix.append(string)

    #some useful properties that can be calculated from the lattice vectors
    def volume(self):
        #calculate the volume of the cell
        return self.a[0]*self.b[1]*self.c[2]-self.a[0]*self.b[2]*self.c[1]+\
                    self.a[1]*self.b[2]*self.c[0]-self.a[1]*self.b[0]*self.c[2]+\
                    self.a[2]*self.b[0]*self.c[1]-self.a[2]*self.b[1]*self.c[0]
    def a_length(self):
        return math.sqrt(self.a[0]**2+self.a[1]**2+self.a[2]**2)
    def b_length(self):
        return math.sqrt(self.b[0]**2+self.b[1]**2+self.b[2]**2)
    def c_length(self):
        return math.sqrt(self.c[0]**2+self.c[1]**2+self.c[2]**2)
    def alpha_angle(self):
        return math.acos(dot(self.b,self.c)/(self.b_length()*self.c_length()))*180/math.pi
    def beta_angle(self):
        return math.acos(dot(self.a,self.c)/(self.a_length()*self.c_length()))*180/math.pi
    def gamma_angle(self):
        return math.acos(dot(self.a,self.b)/(self.a_length()*self.b_length()))*180/math.pi

    #crystalographer's definition, no 2PI
    def reciprocal_a(self):
        V = abs(dot(self.a,cross(self.b,self.c)))
        return [cross(self.b,self.c)[0]/V,cross(self.b,self.c)[1]/V,\
                cross(self.b,self.c)[2]/V]
    def reciprocal_b(self):
        V = abs(dot(self.b,cross(self.c,self.a)))
        return [cross(self.c,self.a)[0]/V,cross(self.c,self.a)[1]/V,\
                cross(self.c,self.a)[2]/V]
    def reciprocal_c(self):
        V = abs(dot(self.c,cross(self.a,self.b)))
        return [cross(self.a,self.b)[0]/V,cross(self.a,self.b)[1]/V,\
                cross(self.a,self.b)[2]/V]

    #physical style defination, 2*PI
    def pi2reciprocal_a(self):
        return [2.0*math.pi*reciprocal_a()[0],2.0*math.pi*reciprocal_a()[1],\
                2.0*math.pi*reciprocal_a()[2]]
    def pi2reciprocal_b(self):
        return [2.0*math.pi*reciprocal_b()[0],2.0*math.pi*reciprocal_b()[1],\
                2.0*math.pi*reciprocal_b()[2]]
    def pi2reciprocal_c(self):
        return [2.0*math.pi*reciprocal_c()[0],2.0*math.pi*reciprocal_c()[1],\
                2.0*math.pi*reciprocal_c()[2]]
    
    def atoms_counts(self):
        if len(self.atoms) != 0:
            kinds = []
            counts = []
            for atom in self.atoms:
                if atom[0] not in kinds:
                    kinds.append(atom[0])
                    counts.append(1)
                else:
                    counts[kinds.index(atom[0])] += 1
        return [kinds,counts]

    def scale_volume(self, scalingfactor):
        #apply scaling to lattice constants
        cuberootscale = scalingfactor**(1.0/3.0)
        for i in range(3):
            self.a[i] *= cuberootscale
            self.b[i] *= cuberootscale
            self.c[i] *= cuberootscale

    def shake(self, dR=0.1):
        #randomly shake the atoms in POSCAR
        for atom in self.atoms:
            for i in range(3):
                atom[1][i] += (-1.0+2.0*random.random())*dR

    #File export options
    def save(self, filename="POSCAR"):
        newfile = file(filename, "w")

        newfile.write(self.description + "\n")

        #put lattice constant to 1.0
        newfile.write("1.00000\n")

        #a,b,c
        a_line = "%2.9f %2.9f %2.9f \n" % (self.a[0],self.a[1],self.a[2])
        b_line = "%2.9f %2.9f %2.9f \n" % (self.b[0],self.b[1],self.b[2])
        c_line = "%2.9f %2.9f %2.9f \n" % (self.c[0],self.c[1],self.c[2])
        newfile.write(a_line)
        newfile.write(b_line)
        newfile.write(c_line)

        #atom counts
        element = map(str, atoms_counts()[0])
        " ".join(element)
        element += "\n"
        newfile.write(element)

        elementnum = map(str, atoms_counts()[1])
        " ".join(element_num)
        element_num += "\n"
        newfile.write(element_num)
        
        #fix information
        if self.constrain:
            newfile.write("Selective\n")

        #coordination types
        if self.direct:
            newfile.write("Direct\n")
        else:
            newfile.write("Cartesian\n")

        #atomic positions
        if self.constrain:
            for i in range(atoms_counts()[1]):
                coord_line = "%2.9f %2.9f %2.9f %s %s %s\n" % (self.atoms[i][1][0],\
                        self.atoms[i][1][1],self.atoms[i][1][2],self.fix[i][0],\
                        self.fix[i][1],self.fix[i][2])
                newfile.write(coord_line)
            newfile.close()
        else:
            for atom in self.atoms:
                coord_line = "%2.9f $2.9f %2.9f \n" % (atom[1][0],atom[1][1],atom[1][2])
                newfile.write(coord_line)
            newfile.close()

class Kpoints():
    def __init__(self, filename="KPOINTS"):
        #read data from KPOINTS file
        kfile = file(filename, 'r')
        kdata = kfile.readlines()
        kfile.close()

        self.gridtype = "unknown"
        self.gridsize = 0
        self.origin = [0.0,0.0,0.0]
        self.description = kdata[0].strip()

        if int(kdata[1][0]) == 0:
            #automatic generation of k-points, as far only support Auto
            if kdata[2][0].upper() == "M":
                #Monkhorst-Pack grid
                self.gridtype = "Monkhorst-Pack"
                self.gridsize = list(map(int, kdata[3].split()))
            elif int(kdata[2][0].upper() == "G"):
                    #Gamma-centred grid
                    self.gridtype = "Gamma"
                    self.gridsize = list(map(int, kdata[3].split()))
        else:
            print("this script only surpport auto-generated k-points now!!!\n")

class Incar():
    #basically a wrapper for the dictionary class        
    def __init__(self, filename="INCAR"):
        #load all input tag into a hash
        self.tags = dict()

        taglines = subprocess.getoutput("grep = "+ filename).split("\n")
        for line in taglines:
            tag_key = line.split("=")[0].strip()
            tag_value = line.split("=")[1].strip()

            #don't load disabled tags
            if tag_key[0] != "#":
                self.tags[tag_key] = tag_value

    def has_tag(self, tagname):
        return tagname in self.tags

    def tag_value(self, tagname):
        return self.tags[tagname]

    def change_tag(self, tagname, value):
        self.tags[tagname] = value
        return tagname in self.tags

    def add_tag(self, tagname, value):
        self.tags[tagname] = value

    def zap_tag(self, tagname):
        self.tags.delete(tagname)

    def save(self, filename="POSCAR"):
        newfile = file(filename, "w")

        taglist = self.tags.keys()
        taglist.sort()

        for tag in taglist:
            newfile.write(tag+" = "+self.tags[tag]+"\n")

        newfile.close()

#some data extract with grep
def get_total_energy(where):
    return float(subprocess.getoutput("grep \"free  energy\" "+where+"|tail -1").split()[4])

def get_ediff(where):
    return float(subprocess.getoutput("grep \"EDIFF \" "+where).split()[2])

def get_scf_delay(where):
    return abs(int(subprocess.getoutput("grep NELMDL "+where).split()[6]))

def get_fermi(where):
    return float(subprocess.getoutput("grep E-fermi "+where).split()[2])

def get_entropy(where):
    return float(subprocess.getoutput("grep EENTRO "+where+"|tail -1").split()[4])

def get_external_pressure(where):
    return float(subprocess.getoutput("grep pressure "+where+"|tail -1").split()[3])

def get_number_of_kpoints(where):
    return int(subprocess.getoutput("grep NKPT "+where+"|tail -1").split()[3])

#Peter's Norving's tiny spell checker in 20 lines... use it to check for misspelled INCAR tags
def words(text):
    return re.findall(r'[a-z]+', text.lower())

def train(features):
    model = collections.defaultdict(lambda: 1)
    for f in features:
        model[f] += 1
    return model

tag_db = train(words("""
  addgrid aexx aggac aggax aldac algo amin amix amix_mag apaco
  bmix bmix_mag cmbj cmbja cmbjb cshift deper dipol dq ebreak
  ediff ediffg efield_pead emax emin enaug encut evenonly
  ferdo ferwe gga gga_compat hfscreen i_constrained_m ialgo
  ibrion icharg ichibare icorelevel images imix inimix iniwav
  ipead isif ismear ispin istart isym iwavr kblock kgamma kpar
  kspacing lambda lasph lcalceps lcalcpol lcharg lchimag lcorr
  ldau ldauj ldaul ldauprint ldautype ldauu ldiag lefg lelf
  lepsilon lhfcalc lhyperfine lkproj lmaxfock lmaxfockae
  lmaxmix lmaxpaw lmaxtau lmixtau lnabla lnmr_sym_red
  lnoncollinear loptics lorbit lpead lplane lreal lrpa
  lscalapack lscalu lsorbit lspectral lthomas lvhar lvtot
  lwannier90 lwannier90_run lwave lwrite_mmn_amn m_constr
  magmom maxmix metagga mixpre nbands nblk nblock ncore ndav
  nedos nelect nelm nelmdl nelmin nfree ngx ngxf ngy ngyf
  ngyromag ngz ngzf nkred nkredx nkredy nkredz nlspline
  nomega nomegar npaco npar nsim nsw nupdown nwrite
  oddonly omegamax omegatl pflat plevel pomass potim
  prec precfock proutine pstress pthreshold quad_efg
  ropt rwigs saxis sigma smass smearings spring symprec
  system tebeg teend time voskown wc weimin zval
  """))

alphabet = 'abcdefghijklmnopqrstuvwxyz'

def edits1(word):
    splits = [(word[:i],word[i:]) for i in range(len(word)+1)]
    deletes = [a + b[1:] for a,b in splits if b]
    transposes = [a + b[1] + b[0] +b[2:] for a,b in splits if len(b)>1]
    replaces = [a + c + b[1:] for a,b in splits for c in alphabet if b]
    inserts = [a + c + b for a,b in splits for c in alphabet]
    return set(deletes + transposes + replaces + inserts)

def known_edits2(word):
    return set(e2 for e1 in edits1(word) for e2 in edits1(e1) if e2 in tag_db)

def known(words):
    return set(w for w in words if w in tag_db)

def correct(word):
    #candidate is set to the first true item
    candidates = known([word]) or known(edits1(word)) or known_edits2(word) or [word]
    return max(candidates, key=tag_db.get)

def read_incar():
    incarfile = open("INCAR","r").readlines()
    result = {}
    for line in incarfile:
        parts = line.strip().split("=")

        if len(parts) >= 2 and parts[0][0] != "#" and parts[0][0] != "!":
            if parts[1].find("#") != -1:
                stripped = parts[1].split("#")[0]
            else:
                stripped = parts[1]
            result[parts[0].strip().upper()] = stripped.strip().upper()
    return result

#count how many magmoms was given on en INCAR lines
def atom_defs_given(data):
    result = 0
    for chunk in data.split():
        if "*" in chunk:
            terms = chunk.split("*")
            if len(terms) != 2:
                logging.warning("Prefight was unable to parse MAGMOM/LDAU, Chunk is %s" % (chunk))
            else:
                result += int(terms[0])
        else:
            result += 1
    return result

#count e.g. how many LDAUs was given on an INCAR line
def atom_kinds_given(data):
    return len(data.split())

parser = OptionParser()
parser.add_option("-o","--output",dest="output",help="Print preflight check data to file",metavar="FILE")
parser.add_option("-l","--loglevel",dest="loglevel",help="Log events of this level and higher")
parser.add_option("-c","--cores",dest="cores",help="How many cores you intend to run on ")
parser.add_option("-n","--nodes",dest="nodes",help="How many nodes you intend to run on")
(options,args) = parser.parse_args()

#TODO: implement LOGLEVEL setting
loglevel = logging.INFO

if options.output:
    logging.basicConfig(filename = options.output,level=loglevel,format='[%(filename)7s] %(message)s')
else:
    logging.basicConfig(level=loglevel,format='[%(levelname)7s] %(message)s')

allfiles = True
#Assert all files exist necessary to run a calculation 
for inputfile in ["INCAR","POSCAR","POTCAR","KPOINTS"]:
    if not os.path.isfile(inputfile):
        allfiles = False
        logging.error("VASP cannot run without the %s file." % (inputfile))

if not allfiles:
    sys.exit(0)

#Read the INCAR and spell check tags
incar = read_incar()
for tag in incar:
    matched = correct(tag.lower())
    if matched != tag.lower():
        logging.warning("Unkoen tag %s in INCAR. Did you mean %s?" % (tag, matched.lower()))

#Look if CONTCAR = POSCAR when restarting
if os.path.isfile("CONTCAR") and subprocess.getoutput("diff CONTCAR POSCAR") != "" and "NSW" in incar:
    logging.warning("POSCAR and CONTCAR are the same. Did you forget to copy CONTCAR to POSCAR?")

#Determine ENMIN and ENMAX from POTCAR
enmax, enmin = 0.0, 0.0
potcarfile = subprocess.getoutput("grep ENMAX POTCAR").split("\n")
for line in potcarfile:
    parts = line.strip().split()
    maxnum = float(parts[2][0:-1])
    minnum = float(parts[5])
    if maxnum > enmax:
        enmax = maxnum
    if minnum > enmin:
        enmin = minnum

#Look if ENCUT is too large or too small
if "ENCUT" not in incar:
    logging.warning("There is no ENCUT value in the INCAR file.")
else:
    cutoff = float(incar["ENCUT"])
    if cutoff < enmin:
        logging.warning("ENCUT is lower than the recommend ENMIN in the POTCAR file.")
    if cutoff < enmax:
        logging.warning("ENCUT is lower than the recommned ENMAX in the POTCAR file. \
                This could lead to artificially compression if you relax the cell?")
    if "ISIF" in incar:
        if (incar["ISIF"] == 3 or incar["ISIF"] == 6 or incar["ISIF"] == 7) and incar["ENCUT"] < enmax*1.3:
            logging.warning("Volume relaxation with ENCUT < ENMAX*1.3. Consider adding PSTRESS.")
    if cutoff > 1.6*enmax:
        logging.warning("ENCUT is mach higher than ENMAX. it might lead to numrical instability.")

#TODO: should abstract acess to tags or check more carefully for existence

#Check for PREC High
if "PREC" in incar and incar["PREC"].lower() == "high":
    logging.warning("The PREC = High setting is deprecated . Use Accurate instead.")

#POTIM should be used with IBRION =1
if "IBRION" in incar and incar["IBRION"] == "1" and "POTIM" not in incar:
    logging.warning("POTIM is not set, even though IBRION=1 is selected. VASP's default value of 0.5 \
            may not be optimal.")

#Check LDA+U and LMAXMIX
if "LDAU" in incar:
    if "LMAXMIX" in incar:
        if int(incar["LMAXMIX"]) < 4:
            logging.info("LDA+U calculated may require higher LMAXMIX if you have d/f-element.")
    else:
        logging.info("LAD+U calcualtions may require higher LMAXMIX if you have d/f-element.")

#Check MAGMOM if spin-polarized
if "ISPIN" in incar and int(incar["ISPIN"]) == 2 and "MAGMOM" not in incar and "ICHARG" in incar \
        and int(incar["ICHARG"]) != 1:
    logging.info("It is a good idea to set MAGMOM for initializing spin-polarized calculations.")

#Check for symmetry restricted MD
if "IBRION" in incar and incar["IBRION"] == 0 and ("ISYM" not in incar or int(incar["ISYM"]) > 0):
    logging.warning("You should turn off symmetry constrains when doing molecular dynamics.")

#parallel relative parameters should be set
if "NPAR" not in incar and "NCORE" not in incar:
    logging.warning("You must set approriate NCORE or NPAR to get good performance.")

#Check number of atoms in cell vs MAGMOM and LDAU parameters
cell = Poscar("POSCAR")
vol = cell.volume()
natoms = sum(cell.atoms_counts()[1])

#Check cell volume, write waning if less than one cubic Angstrom per atom
if vol/natoms < 1.0:
    logging.warning("The volume per atom is less than 1.0 A^3/atom.")

#Check number of elements in MAGMOM and LDA+u, should match atoms.
if "MAGMOM" in incar and atom_defs_given(incar["MAGMOM"]) != natoms:
    logging.error("MAGMOM and number of atopms is inconsistent.")

if "LDAUU" in incar and atom_kinds_given(incar["LDAUU"]) != len(cell.atom_counts()[1]):
    logging.error("LDAUU and number of elements is inconsistent. %d != %d" % (atom_kinds_given(incar["LDAUU"]),\
            len(cell.atom_counts()[1])))

if "LDAUJ" in incar and atom_kinds_given(incar["LDAUJ"]) != len(cell.atom_counts()[1]):
    logging.error("LDAUJ and number of elements is inconsistent. %d != %d" % (atom_kinds_given(incar["LDAUJ"]),\
            len(cell.atom_counts()[1])))

#TODO sugesating optimal parrel paramiters 
