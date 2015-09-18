__author__ = 'Greg Poisson'

import numpy
import math
import math
from sklearn.svm import SVR
from sklearn.linear_model import Lasso
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
    from rdkit.Chem import Fragments
    from rdkit.Chem import AllChem
    from rdkit.Chem import MolSurf

    global featuresFile, numFeatures
    featuresFile = open(fileName, 'w')      # Molecule features output file

    # run gaussian jobs
#    gaussian.setNumMols()
#    gaussian.makeAllGinps()
#    gaussian.runGaussianOnAllGinps()

    # open database file
    drugDB = Chem.SDMolSupplier("FKBP12_binders.sdf")

    if debug:
        print "\n\tNo features data file found. Writing new features data file.\n"

    text = ""       # Placeholder for feature data
    molCount = 0
    convergedCount = 0
    converged_and_different = 0
    drug_name = []

    # load fragment descriptor
    Fragments._LoadPatterns(fileName='/usr/local/anaconda/pkgs/rdkit-2015.03.1-np19py27_0/share/RDKit/Data/FragmentDescriptors.csv')

    # Select features of interest
    for mol in drugDB:
	if molCount > -1:
#		print mol.GetProp("BindingDB Target Chain Sequence")
		gaussian_log_file = "gaussian_files/drug_"+str(molCount)+".log"
		converged, dipole, quadrupole, octapole, hexadecapole, dg_solv = gaussian.parseGaussianLog(gaussian_log_file)
		if converged == "True" and mol.GetProp("BindingDB Target Chain Sequence") == "MGVQVETISPGDGRTFPKRGQTCVVHYTGMLEDGKKFDSSRDRNKPFKFMLGKQEVIRGWEEGVAQMSVGQRAKLTISPDYAYGATGHPGIIPPHATLVFDVELLKLE":
			if convergedCount ==0:
				diff = "True"
			else:
				diff = "True"
				for i in range(converged_and_different):
					if mol.GetProp("BindingDB Ligand Name") == drug_name[i]:
						diff = "False"
						break
			
			if diff == "True":
				drug_name.append(mol.GetProp("BindingDB Ligand Name"))				
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
				text += "{}\n".format(Fragments.fr_Al_OH(mol)) # aliphatic alcohols
				text += "{}\n".format(Fragments.fr_Ar_OH(mol)) # aromatic alcohols
				text += "{}\n".format(Fragments.fr_ketone(mol)) # number of ketones
				text += "{}\n".format(Fragments.fr_ether(mol)) # number of ether oxygens
				text += "{}\n".format(Fragments.fr_ester(mol)) # number of esters
				text += "{}\n".format(Fragments.fr_aldehyde(mol)) # number of aldehydes
				text += "{}\n".format(Fragments.fr_COO(mol)) # number of carboxylic acids
				text += "{}\n".format(Fragments.fr_benzene(mol)) # number of benzenes
		                text += "{}\n".format(Fragments.fr_Ar_N(mol)) # number of aromatic nitrogens
		                text += "{}\n".format(Fragments.fr_NH0(mol)) # number of tertiary amines
		                text += "{}\n".format(Fragments.fr_NH1(mol)) # number of secondary amines
		                text += "{}\n".format(Fragments.fr_NH2(mol)) # number of primary amines
		                text += "{}\n".format(Fragments.fr_amide(mol)) # number of amides
		                text += "{}\n".format(Fragments.fr_SH(mol)) # number of thiol groups
		                text += "{}\n".format(Fragments.fr_nitro(mol)) # number of nitro groups
		                text += "{}\n".format(Fragments.fr_furan(mol)) # number of furan rings
		                text += "{}\n".format(Fragments.fr_imidazole(mol)) # number of imidazole rings
		                text += "{}\n".format(Fragments.fr_oxazole(mol)) # number of oxazole rings
		                text += "{}\n".format(Fragments.fr_morpholine(mol)) # number of morpholine rings
		                text += "{}\n".format(Fragments.fr_halogen(mol)) # number of halogens
				text += "\nKI: {}\n".format(mol.GetProp("Ki (nM)"))
				text += "\n"        # Use a blank line to divide molecule data
				
				featuresFile.write(text)
				text = ""
				converged_and_different += 1
			convergedCount += 1
	else:
		break
       	molCount += 1

    print "Number of molecules with converged gaussian log files and correct sequence:", convergedCount, "\n"
    print "Number of overlap drugs:", convergedCount - converged_and_different
    featuresFile.close()


# Returns 3 arrays; Training data for features, training data for KIs, and all molecules without given KI values
def getFeatures(fileName):
    global featuresFile
    featuresFile = open(fileName, 'r')

    print "Getting features from data file.\n"

    molsWithKI = []
    molsFull = []
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
                if line[1][0] == '>':
                    line[1] = line[1][1:]
                molsWithKI.append(featuresList)
                kiList.append(float(line[1]))
                molsFull.append(featuresList)
            else:
                molsFull.append(featuresList)
        elif f == False:
            f = True

            featuresList = []

    fullKiList = kiList
    training_x = numpy.transpose(numpy.asarray(molsWithKI))
    training_y = numpy.asarray(kiList)
    newData = numpy.transpose(numpy.asarray(molsFull))

    mean = numpy.empty(newData.shape[0],dtype=float)
    std = numpy.empty(newData.shape[0],dtype=float)
    for i in range(newData.shape[0]):
    	mean[i] = numpy.mean(newData[i,:])
    	std[i] = numpy.std(newData[i,:])
        if std[i] == 0:
		newData[i,:] = newData[i,:] - mean[i]
		training_x[i,:] = training_x[i,:] - mean[i]
        else:
        	newData[i,:] = (newData[i,:] - mean[i]) / std[i]
        	training_x[i,:] = (training_x[i,:] - mean[i]) / std[i]


    #print training_x
    #print (training_y)

#    normalized_x = numpy.transpose(normalize(training_x))
#    normalized_y = normalize(training_y)[0]

    #print normalized_x
    #print normalized_y

#    standardized_x = numpy.transpose(scale(training_x))
#    standardized_y = scale(training_y)

    training_x = numpy.transpose(training_x)
    newData = numpy.transpose(newData)

#    assert len(normalized_x) == len(normalized_y)
#    assert len(standardized_x) == len(standardized_y)
#    assert len(normalized_x) == len(training_y)
#    assert len(standardized_x) == len(training_y)

#    if debug:
#        print "{} small molecules in database.".format(len(training_x) + len(newData))
#        print "{} have KI values listed in the database.".format(len(training_y))
#        print "{} do not have KI values listed in the database.\n".format(len(newData))

    return training_x, training_y, newData, numpy.asarray(fullKiList)


def main():

    if debug:
        print "\n\n\tdrugBind.py"

    try:
        train_x, train_y, full_features, full_kis = getFeatures(featuresFilename)
    except IOError:
        makeFeatures(featuresFilename)
        train_x, train_y, full_features, full_kis = getFeatures(featuresFilename)

#    train_y = numpy.log(train_y)

    # machine learning steps
    # fit a SVM model to the data
#    model = SVR(kernel='linear', C=1e10)
    model = SVR(kernel='rbf', C=1e12, gamma=0.1)
    model.fit(train_x, train_y)
    if debug:
        print model
        print "\nUsing training data to test model accuracy:"
    # check
    expected = train_y
    predicted_svr = model.predict(train_x)
    mse = numpy.mean((predicted_svr-expected)**2)
    print("\n\tMean of squared errors for trained SVR: {}".format(mse))
    # make predictions
#    expected = full_kis
#    predicted_svr = model.predict(full_features)

#    print "SVR coefficients:", model.coef_

    # summarize the fit of the model
#    mse = numpy.mean((predicted_svr-expected)**2)
    # mean of squared errors
    if debug:
        print("\n\tMean of squared errors for predicted SVR: {}".format(mse))


    '''
    Returns the coefficient of determination R^2 of the prediction.
    The coefficient R^2 is defined as (1 - u/v), where u is the regression sum of squares
    ((y_true - y_pred) ** 2).sum() and v is the residual sum of squares ((y_true - y_true.mean()) ** 2).sum().
    Best possible score is 1.0, lower values are worse.
    '''
    if debug:
        print("\tModel score: {}".format(model.score(train_x, train_y)))

    print("\nRunning Lasso model...\n")
    model = Lasso()
    model.fit(train_x, train_y)
    print "Lasso coefficients:", model.coef_
    # make predictions
    expected = train_y
    predicted_lasso = model.predict(train_x)
    # summarize the fit of the model
    mse = numpy.mean((predicted_lasso-expected)**2)
    # mean of squared errors
    if debug:
        print("\n\tMean of squared errors for Lasso: {}".format(mse))

    fit_out = open("actual_v_predicted_kis.rbf.KI_lt_2_model.full_data.dat", "w")
    for i in range(expected.size):
        fit_out.write("%20.12f %20.12f %20.12f\n" % (predicted_svr[i], predicted_lasso[i], expected[i]))
    fit_out.close()

main()
