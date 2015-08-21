__author__ = 'Greg Poisson'

import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

ginps = []
numMols = 0
drugDB = None

def setNumMols(dbName="FKBP12_binders.sdf"):
    global numMols, drugDB
    drugDB = Chem.SDMolSupplier(dbName)
    lastMol = None

    for mol in drugDB:
        if mol != lastMol:
            lastMol = mol
            ginpName = "drug_{}.ginp".format(numMols)
            ginps.append(ginpName)
            numMols += 1

def makeAllGinps(dbName="FKBP12_binders.sdf"):
    global ginps

    count = 0
    lastMol = None
    for mol in drugDB:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        block = Chem.MolToMolBlock(mol)     # generates string block of molecule coord/bond data
        if mol != lastMol:
            lastMol = mol
            ginpName = "drug_{}.ginp".format(count)
            ginp = open(ginpName, 'w')
            text = "%nprocshared=3\n"
            text += "%mem=4GB\n"
            text += "%Chk=drug_{}.chk\n".format(count)
            text += "# B3LYP/6-31G(d,p) opt pop=(MK,saveESP)\n\n"
            text += "Optimize and compute dipole, quadrupole, etc...\n\n"

            text += " 0  1\n"

            bCount = 0
            line = ""
            for ch in block:
                if ch == "\n":
                    data = line.split()
                    if (bCount >= 4) & (len(data) == 16):
                        element = data[3]
                        coords = "\t {0:.6f}\t{0:.6f}\t{0:.6f}\n".format(float(data[0]), float(data[1]), float(data[2]))
                        text += element + coords
                    line = ""
                else:
                    line += ch
                bCount += 1

            text += "\n--Link1--\n"
            text += "%Chk=drug {}.chk\n".format(count)
            text += "%NoSave\n"
            text += "# B3LYP/6-31G(d) Geom=Check Guess=Read SCRF=PCM\n\n"
            text += "Compute solvated energy\n"
            text += " 0 1\n"

            ginp.write(text)
            ginp.close()

            count += 1

            sys.stdout.write("Writing Gaussian input files: {0:.2f}% Complete\r".format((float(count) / float(numMols)) * 100))
            sys.stdout.flush()


def runGaussianOnAllGinps():
    count = 0
    for ginp in ginps:
        if count < 5:
            os.system("g09 {} > {}.log".format(ginp, ginp[:-5]))
            count += 1
            sys.stdout.write("Running Gaussian on all inputs: {0:.2f}% Complete\r".format((float(count) / float(numMols)) * 100))
            sys.stdout.flush()

