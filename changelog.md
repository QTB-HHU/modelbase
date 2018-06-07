# Changelog
## 0.2.1
* Changed total concentration names in label models to base_name+"_total"
* Added print_stoichiometryMatrix function to return pandas dataframe
* Added set_reaction and set_reaction_v shortcut functions
* Allowed ParameterSet input as update method for parameter update
* Warning when overwriting parameters with ParameterSet.update()  


## 0.2.0
* Separated analysis methods from model class to analysis.py
* Unified Simulator class calls with constructor method modelbase.Simulator()
* Removed AlgmSimulate class
* Unified AlgmModel and Model classes
* Changed algebraic module construction


## 0.1.8
* bugfix: in LabelModel setting c=0 (no labels in this compound) led to an error, because the sum of all labels
had the same name as the compound. Fixed.
* verbosity can be passed to the assimulo solver.


## 0.1.7
Support for the differential equation solver sundials (CVODE)
through the python package [assimulo](http://www.jmodelica.org/assimulo).
Falls back to scipy.integrate.ode if assimulo cannot be loaded.

Brief installation instructions of sundials/assimulo (tested on Ubuntu 14.04 and 16.04 and MacOX X):
* Install sundials-2.6.0 [here](https://computation.llnl.gov/projects/sundials/sundials-software). The version is important. We did not get 3.0.0 to run. You will need cmake and ccmake for it.
Set compiler flag -fPIC.
* pip install assimulo