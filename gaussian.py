__author__ = 'Greg Poisson'

import os
import sys
import os.path
import numpy
import math
from rdkit import Chem
from rdkit.Chem import AllChem

ginps = []
numMols = 0
drugDB = None

def setNumMols(dbName="FKBP12_binders.sdf"):
    global numMols, drugDB
    drugDB = Chem.SDMolSupplier(dbName)
    lastMol = None

    count=0
    for mol in drugDB:
        if mol != lastMol:
            lastMol = mol
            print "Molecule: ", count+1, "mol.GetProp(Ki (nM))", mol.GetProp("Ki (nM)"), "\n"
	    try:
		mol.GetProp("Ki (nM)")
            except KeyError, e:	
	        pass
            else:
	        print "Molecule ", count+1, "has a KI\n"
                ginpName = "drug_{}.ginp".format(numMols)
                ginps.append(ginpName)
                numMols += 1
            count += 1

    print "Number of Molecules:", numMols, "\n"

def makeAllGinps(dbName="FKBP12_binders.sdf"):
    global ginps

    count = 0
    lastMol = None
    for mol in drugDB:
        if (count > 0 and count < 500):
            mol = Chem.AddHs(mol)
            if mol != lastMol:
                lastMol = mol
                try:
            	    mol.GetProp("Ki (nM)")
                except KeyError,e:
                    pass
                else:
                    ginpName = "gaussian_files/drug_{}.ginp".format(count)
                    chkName = "gaussian_files/drug_{}.chk".format(count)
#                    if os.path.isfile(ginpName) == "False":
		    if count > 0:
            	        AllChem.EmbedMolecule(mol)
            	        AllChem.UFFOptimizeMolecule(mol)
            	        block = Chem.MolToMolBlock(mol)     # generates string block of molecule coord/bond data
                        ginp = open(ginpName, 'w')
                        text = "%nprocshared=5\n"
                        text += "%mem=4GB\n"
                        text += "%Chk=drug_{}.chk\n".format(count)
			if os.path.isfile(chkName) == "False":
                        	text += "# HF/6-31G opt=(maxCycles=50)\n\n"
			else:
                        	text += "# HF/6-31G geom=Check Guess=Read opt=(maxCycles=50)\n\n"
                        text += "Optimize at HF/6-31G level of theory\n\n"
     
                        text += " 0  1\n"
     
			if os.path.isfile(chkName) == "False":
                            bCount = 0
                            line = ""
                            for ch in block:
                                if ch == "\n":
                                    data = line.split()
                                    if (bCount >= 4) & (len(data) == 16):
                                        element = data[3]
                                        #coords = "\t {0:.6f}\t{0:.6f}\t{0:.6f}\n".format(float(data[0]), float(data[1]), float(data[2]))
                                        coords = " " + str(data[0]) + " " + str(data[1]) + " " + str(data[2]) + "\n"
                                        text += element + coords
                                    line = ""
                                else:
                                    line += ch
                                bCount += 1
     
                        text += "\n--Link1--\n"
                        text += "%nprocshared=5\n"
                        text += "%mem=4GB\n"
                        text += "%Chk=drug_{}.chk\n".format(count)
                        text += "# B3LYP/6-31G(d,p) opt Geom=Check pop=(MK,saveESP)\n\n"
                        text += "Compute solvated energy\n\n"
                        text += " 0  1\n"
                        text += "\n--Link1--\n"
                        text += "%nprocshared=5\n"
                        text += "%mem=4GB\n"
                        text += "%Chk=drug_{}.chk\n".format(count)
                        text += "%NoSave\n"
                        text += "# B3LYP/6-31G(d,p) Geom=Check Guess=Read SCRF=PCM\n\n"
                        text += "Compute solvated energy\n\n"
                        text += " 0  1\n\n"
     
                        ginp.write(text)
                        ginp.close()

        count += 1

        sys.stdout.write("Writing Gaussian input files: {0:.2f}% Complete\r".format((float(count) / float(numMols)) * 100))
        sys.stdout.flush()


def runGaussianOnAllGinps():
    count = 0
    for ginp in ginps:
        if (count > 0 and count < 500):
#         print "Optimizing drug ", count, "\n"
             os.system("cd gaussian_files; g09 {} > {}.log; rm -f *.rwf ;cd ..".format(ginp, ginp[:-5]))
             sys.stdout.write("Running Gaussian on all inputs: {0:.2f}% Complete\r".format((float(count) / float(500)) * 100))
             sys.stdout.flush()
        count += 1

def parseGaussianLog(fileName):

	glog = 	open(fileName,"r")
	glog_lines = glog.readlines()
	glog.close()

	converged = checkConvergence(glog_lines)

	if converged == "True":
		print "Converged\n"
		dipole = getDipole(glog_lines)
		quadrupole = getQuadrupole(glog_lines)
		octapole = getOctapole(glog_lines)
		hexadecapole = getHexadecapole(glog_lines)
		DGSolv = getDGSolv(glog_lines)
	else:
		print "Did not converge\n"
                dipole = 0
                quadrupole = 0
                octapole = 0
                hexadecapole = 0
                DGSolv = 0


	return dipole, quadrupole, octapole, hexadecapole, DGSolv

# determine if geometry optimizations in gaussian converged.  Note this checks for two subsequent geometry optimizations
def checkConvergence(log_lines):

	first_convergence_flag = "False"
	second_convergence_flag = "False"

	for i in range(len(log_lines)):
		if "        Item               Value     Threshold  Converged?" in log_lines[i]:
			next_line = log_lines[i+1]
			array = next_line.split()
			max_force_convergence = array[4]
			next_line = log_lines[i+2]
			array = next_line.split()
			rms_force_convergence = array[4]
			next_line = log_lines[i+3]
			array = next_line.split()
			max_disp_convergence = array[4]
			next_line = log_lines[i+4]
			array = next_line.split()
			rms_disp_convergence = array[4]
			if first_convergence_flag == "False" and max_force_convergence == "YES" and rms_force_convergence == "YES" and max_disp_convergence == "YES" and rms_disp_convergence == "YES":
				first_convergence_flag = "True"
			elif max_force_convergence == "YES" and rms_force_convergence == "YES" and max_disp_convergence == "YES" and rms_disp_convergence == "YES":
				second_convergence_flag = "True"
		if second_convergence_flag == "True":
			break

	if first_convergence_flag == "True" and second_convergence_flag == "True":
		return "True"
	else:
		return "False"


#routine to parse a gaussian log file and extract the magnitude of the last dipole moment read
def getDipole(log_lines):

	dipoleFlag = "False"
	for i in range(len(log_lines)-1,-1,-1):
		if "Dipole moment (field-independent basis, Debye):" in log_lines[i]:
			dipoleFlag = "True"
			next_line = log_lines[i+1]
			array = next_line.split()
			dipole = float(array[7])
			break

	if dipoleFlag == "True":
		pass
	else:
		dipole = 0

	return dipole


#routine to parse a gaussian log file and extract the magnitude of the last quadrupole moment read
def getQuadrupole(log_lines):

	quadrupoleFlag = "False"
	quadrupole = 0
	for i in range(len(log_lines)-1,-1,-1):
		if "Quadrupole moment (field-independent basis, Debye-Ang):" in log_lines[i]:
			quadrupoleFlag = "True"
			next_line = log_lines[i+1]
			array = next_line.split()
			quadrupole += float(array[1])*float(array[1]) + float(array[3])*float(array[3]) + float(array[5])*float(array[5])
                        next_line = log_lines[i+2]
                        array = next_line.split()
                        quadrupole += float(array[1])*float(array[1]) + float(array[3])*float(array[3]) + float(array[5])*float(array[5])
			quadrupole = math.sqrt(quadrupole)
			break

	if quadrupoleFlag == "True":
		pass
	else:
		quadrupole = 0

	return quadrupole

#routine to parse a gaussian log file and extract the magnitude of the last octapole moment read
def getOctapole(log_lines):

	octapoleFlag = "False"
	octapole = 0
	for i in range(len(log_lines)-1,-1,-1):
		if "Octapole moment (field-independent basis, Debye-Ang**2):" in log_lines[i]:
			octapoleFlag = "True"
			next_line = log_lines[i+1]
			array = next_line.split()
			octapole += float(array[1])*float(array[1]) + float(array[3])*float(array[3]) + float(array[5])*float(array[5]) + float(array[7])*float(array[7])
                        next_line = log_lines[i+2]
                        array = next_line.split()
			octapole += float(array[1])*float(array[1]) + float(array[3])*float(array[3]) + float(array[5])*float(array[5]) + float(array[7])*float(array[7])
                        next_line = log_lines[i+3]
                        array = next_line.split()
                        octapole += float(array[1])*float(array[1]) + float(array[3])*float(array[3])
			octapole = math.sqrt(octapole)
			break

	if octapoleFlag == "True":
		pass
	else:
		octapole = 0

	return octapole

#routine to parse a gaussian log file and extract the magnitude of the last octapole moment read
def getHexadecapole(log_lines):

	hexadecapoleFlag = "False"
	hexadecapole = 0
	for i in range(len(log_lines)-1,-1,-1):
		if "Hexadecapole moment (field-independent basis, Debye-Ang**3):" in log_lines[i]:
			hexadecapoleFlag = "True"
			next_line = log_lines[i+1]
			array = next_line.split()
			hexadecapole += float(array[1])*float(array[1]) + float(array[3])*float(array[3]) + float(array[5])*float(array[5]) + float(array[7])*float(array[7])
                        next_line = log_lines[i+2]
                        array = next_line.split()
			hexadecapole += float(array[1])*float(array[1]) + float(array[3])*float(array[3]) + float(array[5])*float(array[5]) + float(array[7])*float(array[7])
                        next_line = log_lines[i+3]
                        array = next_line.split()
			hexadecapole += float(array[1])*float(array[1]) + float(array[3])*float(array[3]) + float(array[5])*float(array[5]) + float(array[7])*float(array[7])
                        next_line = log_lines[i+4]
                        array = next_line.split()
                        hexadecapole += float(array[1])*float(array[1]) + float(array[3])*float(array[3]) + float(array[5])*float(array[5])
			hexadecapole = math.sqrt(hexadecapole)
			break

	if hexadecapoleFlag == "True":
		pass
	else:
		hexadecapole = 0

	return hexadecapole


#routine to parse a gaussian log file and extract the DU of solvation
def getDGSolv(log_lines):

	count = 0
	energy =  numpy.empty(2,dtype=float)
	for i in range(len(log_lines)-1,-1,-1):
		if "SCF Done:  E(RB3LYP) =" in log_lines[i]:
			array = log_lines[i].split()
			energy[count] = float(array[4])
			count += 1
			if count > 1:
				break


	return (energy[0]-energy[1])*627.5095


