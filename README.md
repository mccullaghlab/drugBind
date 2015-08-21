# drugBind
Machine learning python code for predicting binding affinity

If you do not have a features.dat file (there is one posted alongside this code), put an SDF database file in the program directory and drugBind.py will compile a features.dat file for you. If you wish to change the feature set, make the appropriate changes to "makeFeatures()", delete the current features.dat file and run drugBind.py again with the small molecule SDF database in your program directory.

With your features.dat file in the program directory, drugBind.py will import the features into numpy arrays for use in scikit learn using a support vector machine for linear regression analysis.

Use:

    ~$: python drugBind.py
