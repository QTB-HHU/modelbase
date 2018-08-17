# modelbase

modelbase is a python package to help you build and analyze dynamic mathematical models of biological systems. It has originally been designed for the simulation of metabolic systems, but can be used for virtually any processes, in which some substances get converted into others.

modelbase incorporates an easy construction method to define 'reactions'. A rate law and the stoichiometry need to be specified, and the system of differential equations is assembled automatically.

modelbase allows 'algebraic modules', which are useful to implement rapid equilibrium or quasi steady-state approximations. In the simplest instance, they allow easy incorporation of conserved quantities.

modelbase also allows a simple construction of isotope-specific models. This class contains a constructor method that automatically construct all isotope specific versions of a particular reaction. Very cool - check it out!

## Release notes
Version 0.2.3 is the official release for the submission of the
mansucript "Building mathematical models of biological systems
with modelbase, a Python package for semi-automatic ODE assembly
and construction of isotope-specific models" to the Journal of Open
Research Software.

See changelog.md for details on changes of earlier versions.

## Dependencies
* NumPy
* Scipy
* Matplotlib
* Assimulo
* Sundials

## Installation

```python3
pip install modelbase
```

To enable assimulo support, either conda install the package

```python3
conda install -c conda-forge assimulo
# or
conda install -c chria assimulo
```

or build sundials from source using our [installation guide](https://gitlab.com/ebenhoeh/modelbase/blob/master/docs/sundials-installation.md).


## License
[GPL 3](https://gitlab.com/ebenhoeh/modelbase/blob/master/LICENSE)

## Documentation
The official documentation is hosted on [ReadTheDocs](https://modelbase.readthedocs.io/en/latest/).

## Issues and support
If you experience issues using the software please contact us through our [issues](https://gitlab.com/ebenhoeh/modelbase/issues) page. For usage questions, feel free to contact @ebenhoeh .

## Contributing to modelbase
All contributions, bug reports, bug fixes, documentation improvements, enhancements and ideas are welcome.
Currently there is no established contribution mechanism, but any help fixing the issues on our [issues](https://gitlab.com/ebenhoeh/modelbase/issues) page is appreciated. If you found a solution for any problem, feel free to contact @ebenhoeh .
