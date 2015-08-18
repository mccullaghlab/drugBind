# drugBind
Machine learning python code for predicting binding affinity

If you do not have a features.dat file (there is one posted alongside this code), put an SDF database file in the program directory and drugBind.py will compile a features.dat file for you.

With your features.dat file in the program directory, drugBind.py will import the features into numpy arrays for use in scikit learn using a support vector machine for linear regression analysis.

Use:

    ~$: python drugBind.py
