import openbabel
import pybel
import sys 
import re
import os

import itertools 

#####################################################################

def add_F (mol, aidx):

  newid = mol.NumAtoms()
  aatom = openbabel.OBAtom()
  aatom.SetAtomicNum(9)
  mol.AddAtom(aatom)
  mol.AddBond(aidx, newid+1, 1)

  return

#####################################################################

def add_H (mol, aidx):

  newid = mol.NumAtoms()
  aatom = openbabel.OBAtom()
  aatom.SetAtomicNum(1)
  mol.AddAtom(aatom)
  mol.AddBond(aidx, newid+1, 1)

  return

#####################################################################

def add_Cl (mol, aidx):

  newid = mol.NumAtoms()
  aatom = openbabel.OBAtom()
  aatom.SetAtomicNum(17)
  mol.AddAtom(aatom)
  mol.AddBond(aidx, newid+1, 1)

  return

#####################################################################

def add_I (mol, aidx):

  newid = mol.NumAtoms()
  aatom = openbabel.OBAtom()
  aatom.SetAtomicNum(53)
  mol.AddAtom(aatom)
  mol.AddBond(aidx, newid+1, 1)

  return

#####################################################################

def add_Br (mol, aidx):

  newid = mol.NumAtoms()
  aatom = openbabel.OBAtom()
  aatom.SetAtomicNum(35)
  mol.AddAtom(aatom)
  mol.AddBond(aidx, newid+1, 1)

  return

#####################################################################

def add_OH (mol, aidx):

  newid = mol.NumAtoms()

  aatom = openbabel.OBAtom()
  aatom.SetAtomicNum(8)
  mol.AddAtom(aatom)
  mol.AddBond(aidx, newid+1, 1)

  aatom1 = openbabel.OBAtom()
  aatom1.SetAtomicNum(1)
  mol.AddAtom(aatom1)
  mol.AddBond(newid+1, newid+2, 1)

  return

#####################################################################

def add_SH (mol, aidx):

  newid = mol.NumAtoms()

  aatom = openbabel.OBAtom()
  aatom.SetAtomicNum(16)
  mol.AddAtom(aatom)
  mol.AddBond(aidx, newid+1, 1)

  aatom1 = openbabel.OBAtom()
  aatom1.SetAtomicNum(1)
  mol.AddAtom(aatom1)
  mol.AddBond(newid+1, newid+2, 1)

  return

#####################################################################

def add_CH3 (mol, aidx):

  newid = mol.NumAtoms()

  aatom = openbabel.OBAtom()
  aatom.SetAtomicNum(6)
  mol.AddAtom(aatom)
  mol.AddBond(aidx, newid+1, 1)

  aatom1 = openbabel.OBAtom()
  aatom1.SetAtomicNum(1)
  mol.AddAtom(aatom1)
  mol.AddBond(newid+1, newid+2, 1)

  aatom2 = openbabel.OBAtom()
  aatom2.SetAtomicNum(1)
  mol.AddAtom(aatom2)
  mol.AddBond(newid+1, newid+3, 1)

  aatom3 = openbabel.OBAtom()
  aatom3.SetAtomicNum(1)
  mol.AddAtom(aatom3)
  mol.AddBond(newid+1, newid+4, 1)

  return

#####################################################################

def add_CF3 (mol, aidx):

  newid = mol.NumAtoms()

  aatom = openbabel.OBAtom()
  aatom.SetAtomicNum(6)
  mol.AddAtom(aatom)
  mol.AddBond(aidx, newid+1, 1)

  aatom1 = openbabel.OBAtom()
  aatom1.SetAtomicNum(9)
  mol.AddAtom(aatom1)
  mol.AddBond(newid+1, newid+2, 1)

  aatom2 = openbabel.OBAtom()
  aatom2.SetAtomicNum(9)
  mol.AddAtom(aatom2)
  mol.AddBond(newid+1, newid+3, 1)

  aatom3 = openbabel.OBAtom()
  aatom3.SetAtomicNum(9)
  mol.AddAtom(aatom3)
  mol.AddBond(newid+1, newid+4, 1)

  return

#####################################################################

def add_NH2 (mol, aidx):

  newid = mol.NumAtoms()

  aatom = openbabel.OBAtom()
  aatom.SetAtomicNum(7)
  mol.AddAtom(aatom)
  mol.AddBond(aidx, newid+1, 1)

  aatom1 = openbabel.OBAtom()
  aatom1.SetAtomicNum(1)
  mol.AddAtom(aatom1)
  mol.AddBond(newid+1, newid+2, 1)

  aatom2 = openbabel.OBAtom()
  aatom2.SetAtomicNum(1)
  mol.AddAtom(aatom2)
  mol.AddBond(newid+1, newid+3, 1)

  return

#####################################################################

def add_OCH3 (mol, aidx):

  newid = mol.NumAtoms()

  aatom = openbabel.OBAtom()
  aatom.SetAtomicNum(8)
  mol.AddAtom(aatom)
  mol.AddBond(aidx, newid+1, 1)

  aatom1 = openbabel.OBAtom()
  aatom1.SetAtomicNum(6)
  mol.AddAtom(aatom1)
  mol.AddBond(newid+1, newid+2, 1)

  aatom3 = openbabel.OBAtom()
  aatom3.SetAtomicNum(1)
  mol.AddAtom(aatom3)
  mol.AddBond(newid+2, newid+3, 1)

  aatom4 = openbabel.OBAtom()
  aatom4.SetAtomicNum(1)
  mol.AddAtom(aatom4)
  mol.AddBond(newid+2, newid+4, 1)

  aatom5 = openbabel.OBAtom()
  aatom5.SetAtomicNum(1)
  mol.AddAtom(aatom5)
  mol.AddBond(newid+2, newid+5, 1)

  return

#####################################################################

def add_NO2 (mol, aidx):

  newid = mol.NumAtoms()

  aatom = openbabel.OBAtom()
  aatom.SetAtomicNum(7)
  mol.AddAtom(aatom)
  mol.AddBond(aidx, newid+1, 1)

  aatom1 = openbabel.OBAtom()
  aatom1.SetAtomicNum(8)
  mol.AddAtom(aatom1)
  mol.AddBond(newid+1, newid+2, 2)

  aatom3 = openbabel.OBAtom()
  aatom3.SetAtomicNum(8)
  mol.AddAtom(aatom3)
  mol.AddBond(newid+1, newid+3, 1)

  aatom4 = openbabel.OBAtom()
  aatom4.SetAtomicNum(1)
  mol.AddAtom(aatom4)
  mol.AddBond(newid+3, newid+4, 1)

  return

#####################################################################

def mk3dandwrite (mol, i, t3d):

  if (t3d):
    mol.make3D(forcefield='mmff94', steps=50)

  #mol.draw(show=False, update=True)

  outfname = str(i)+".sdf"
  if os.path.exists(outfname):
    os.remove(outfname)
  
  mol.write('sdf', outfname)

#####################################################################

if (len(sys.argv)) == 3:
  filename = sys.argv[1]
else:
  print "usage :", sys.argv[0] , "\"n1;n2;...\" molfile"
  exit(1)

filename = sys.argv[2]
atomnumlist = sys.argv[1]

p = re.compile(r'\s+')
atomnumlist = p.sub(' ', atomnumlist)
atomnumlist = atomnumlist.lstrip()
atomnumlist = atomnumlist.rstrip()

alist = atomnumlist.split(";")

atmset = []

for molin in pybel.readfile( "sdf", filename):
  print "AtomList: "
  for a in alist:
    atomnum = int(a)

    if atomnum < molin.OBMol.NumAtoms():
      print "  " , atomnum
      atmset.append(atomnum)

add_res = {"A"  : add_H,
           "B"  : add_F,
           "C"  : add_Cl,
           "D"  : add_I,
           "E"  : add_Br,
           "F"  : add_OH,
           "G"  : add_CH3,
           "H"  : add_NH2,
           "I"  : add_CF3,
           "L"  : add_SH,
           "M"  : add_OCH3,
           "N"  : add_NO2
           }

outfname = "notunique.sdf"
if os.path.exists(outfname):
  os.remove(outfname)

outfname = "unique.sdf"
if os.path.exists(outfname):
  os.remove(outfname)

inchis = []
molid = 0
for molin in pybel.readfile( "sdf", filename):
  for x in itertools.product('ABCDEFGHILMN', repeat=len(atmset)):
    combo = {}
    for i in range(len(x)):
      combo[atmset[i]] = x[i]

    print molid, combo

    clone = pybel.ob.OBMol(molin.OBMol)
    mol = pybel.Molecule(clone)

    for aid, res in combo.iteritems(): 
      add_res[res](mol.OBMol, aid)

    mol.addh()
    inchicode = mol.write("inchi")
    outfname = str(molid)+".sdf"

    if not (inchicode in inchis):
      inchis.append(inchicode)
      mk3dandwrite(mol, molid, False)

      with open("unique.sdf", "a") as myfile:
        for line in open(str(str(molid)+".sdf")):
           myfile.write(line)

    mk3dandwrite(mol, molid, False)
    with open("notunique.sdf", "a") as myfile:
      for line in open(str(str(molid)+".sdf")):
         myfile.write(line)

    if os.path.exists(outfname):
      os.remove(outfname)

    molid = molid + 1
