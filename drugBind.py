__author__ = 'Greg Poisson'

import numpy
from sklearn.svm import SVR
from sklearn.preprocessing import normalize
from sklearn.preprocessing import scale
import gaussian


debug = True
featuresFilename = "features.dat"
KIFilename = "ki.dat"
featuresFile = None
kiFile = None
features = None
ki = None

# Builds features file from sdf database
def makeFeatures(fileName):


    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import MolSurf

    global featuresFile, numFeatures
    featuresFile = open(fileName, 'w')      # Molecule features output file

    # run gaussian jobs
    gaussian.setNumMols()
    gaussian.makeAllGinps()
    gaussian.runGaussianOnAllGinps()

    # open database file
    drugDB = Chem.SDMolSupplier("FKBP12_binders.sdf")

    if debug:
        print "\n\tNo features data file found. Writing new features data file.\n"

    text = ""       # Placeholder for feature data
    molCount = 0

    # Select features of interest
    for mol in drugDB:
	gaussian_log_file = "gaussian_files/drug_"+str(molCount)+".log"
	dipole, quadrupole, octapole, hexadecapole, dg_solv = gaussian.parseGaussianLog(gaussian_log_file)
        #text += "{}\n".format(molCount)
        text += "{}\n".format(AllChem.ComputeMolVolume(mol))
        text += "{}\n".format(MolSurf.pyLabuteASA(mol))
        text += "{}\n".format(mol.GetNumAtoms())
        text += "{}\n".format(mol.GetNumBonds())
        text += "{}\n".format(mol.GetNumHeavyAtoms())
	text += "{}\n".format(dipole)
	text += "{}\n".format(quadrupole)
	text += "{}\n".format(octapole)
	text += "{}\n".format(hexadecapole)
	text += "{}\n".format(dg_solv)
        text += "\nKI: {}\n".format(mol.GetProp("Ki (nM)"))

        text += "\n"        # Use a blank line to divide molecule data

        featuresFile.write(text)
        text = ""

        molCount += 1

    featuresFile.close()


# Returns 3 arrays; Training data for features, training data for KIs, and all molecules without given KI values
def getFeatures(fileName):
    global featuresFile
    featuresFile = open(fileName, 'r')

    print "Getting features from data file.\n"

    molsWithKI = []
    molsWithoutKI = []
    featuresList = []
    kiList = []

    f = True

    for line in featuresFile:
        line = line.split()
        if f & (len(line) >= 1):    # Reading a valid feature
            featuresList.append(float(line[0]))
        elif f & (len(line) < 1):
            f = False
        elif (f == False) & (len(line) >= 1):
            if len(line) > 1:
                molsWithKI.append(featuresList)
                if line[1][0] == '>':
                    line[1] = line[1][1:]
                kiList.append(float(line[1]))
            else:
                molsWithoutKI.append(featuresList)
        elif f == False:
            f = True

            featuresList = []

    training_x = numpy.transpose(numpy.asarray(molsWithKI))
    training_y = numpy.asarray(kiList)

    #print training_x
    #print (training_y)

    normalized_x = numpy.transpose(normalize(training_x))
    normalized_y = normalize(training_y)[0]

    #print normalized_x
    #print normalized_y

    standardized_x = numpy.transpose(scale(training_x))
    standardized_y = scale(training_y)

    training_x = numpy.transpose(training_x)

    newData = numpy.asarray(molsWithoutKI)

    assert len(normalized_x) == len(normalized_y)
    assert len(standardized_x) == len(standardized_y)

    if debug:
        print "{} small molecules in database.".format(len(training_x) + len(newData))
        print "{} have KI values listed in the database.".format(len(normalized_y))
        print "{} do not have KI values listed in the database.\n".format(len(newData))

    return normalized_x, training_y, newData


def main():

    if debug:
        print "\n\n\tdrugBind.py"

    try:
        train_x, train_y, newData = getFeatures(featuresFilename)
    except IOError:
        makeFeatures(featuresFilename)
        train_x, train_y, newData = getFeatures(featuresFilename)

    # machine learning steps
    # fit a SVM model to the data
    model = SVR()
    model.fit(train_x, train_y)
    if debug:
        print model
        print "\nUsing training data to test model accuracy:"

    # make predictions
    expected = train_y
    predicted = model.predict(train_x)

    # summarize the fit of the model
    mse = numpy.mean((predicted-expected)**2)
    # mean of squared errors
    if debug:
        print("\n\tMean of squared errors: {}".format(mse))


    '''
    Returns the coefficient of determination R^2 of the prediction.
    The coefficient R^2 is defined as (1 - u/v), where u is the regression sum of squares
    ((y_true - y_pred) ** 2).sum() and v is the residual sum of squares ((y_true - y_true.mean()) ** 2).sum().
    Best possible score is 1.0, lower values are worse.
    '''
    if debug:
        print("\tModel score: {}".format(model.score(train_x, train_y)))

main()
