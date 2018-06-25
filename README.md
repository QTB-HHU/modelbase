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

Version 0.2.3 is the official release for the submission of the 
mansucript "Building mathematical models of biological systems 
with modelbase, a Python package for semi-automatic ODE assembly 
and construction of isotope-specific models" to the Journal of Open 
Research Software.

See changelog.md for details on changes of earlier versions.
