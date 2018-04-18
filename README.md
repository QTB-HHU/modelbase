# Modelbase

Modelbase is a python package to help you build and analyse dynamic mathematical models
of biological systems.
It has originally been designed for the simulation of metabolic systems, but can be used
for virtually any processes, in which some substances get converted into others.

Modelbase incorporates an easy construction method to define 'reactions'. A rate law
and the stoichiometry need to be specified, and the system of differential equations
is assembled automatically.

Modelbase allows 'algebraic modules', which are useful to implement rapid equilibrium
or quasi steady-state approximations. In the simplest instance, they allow easy 
incorporation of conserved quantities.

Modelbase also allows a simple construction of isotope-specific models. This class
contains a constructor method that automatically construct all iosotope specific versions
of a particular reaction. Very cool - check it out!

## Installation

```
pip install modelbase
```

## Release notes

Version 0.1.8

Two minor changes:
1. bugfix: in LabelModel setting c=0 (no labels in this compound) led to an error, because the sum of all labels
had the same name as the compound. Fixed.
2. verbosity can be passed to the assimulo solver.


Version 0.1.7

Support for the differential equation solver sundials (CVODE)
through the python package [assimulo](http://www.jmodelica.org/assimulo).
Falls back to scipy.integrate.ode if assimulo cannot be loaded.

Brief installation instructions of sundials/assimulo (tested on Ubuntu 14.04 and 16.04 and MacOX X):
* Install sundials-2.6.0 [here](https://computation.llnl.gov/projects/sundials/sundials-software). The version is important. We did not get 3.0.0 to run. You will need cmake and ccmake for it.
Set compiler flag -fPIC.
* pip install assimulo


