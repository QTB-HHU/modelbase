# This file is part of modelbase.
#
# modelbase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# modelbase is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with modelbase.  If not, see <http://www.gnu.org/licenses/>.

"""Model

Description of the module

"""
import numpy as np
import itertools
import re
import pickle
import pprint
import pandas as pd

from collections import defaultdict
from .parameters import ParameterSet

class Model(object):
    """The main class for modeling. Provides model construction tools.

    Attributes
    ----------
    par : modelbase.parameters.ParameterSet
        Model parameters, stored in a ParameterSet object
    cpdNames : list
        List of all compounds in the model
    cpdIdDict : dict
        Dictionary with all compounds and corresponding IDs
    rateFn : dict
        Contains all reaction names as keys and reaction functions as values
    stoichiometries : dict
        Contains all reaction names as keys and stoichiometry dictionaries as values
    algebraicModules : dict
        Contains all algebraicModule names as keys and algebraic functions as values

    Methods
    -------
    rateNames()
        Returns a list of the rate names
    set_cpds(cpdList)
        Sets the compounds to the given list of compounds
    add_cpd(cpdName)
        Adds a compound to the model
    add_cpds(cpdList)
        Adds a list of compounds to the model
    cpdIds()
        Returns a dictionary with the compound names as keys and their index as values
    find_re_argids(regexp)
        Returns indices for compound names matching the regular expression
    set_rate(rateName, fn, *args)
        Sets a rate
    set_ratev(rateName, fn, *args)
        Sets a rate that depends on additional information (time)
    set_stoichiometry(rateName, stDict)
        Sets the reaction stoichiometry for a rate
    set_stoichiometry_byCpd
        Sets the reaction stoichiometry of a rate by compound name
    add_reaction(rateName, fn, stDict, *args)
        Adds a reaction containing a rate function and stoichiometry
    add_reaction_v(rateName, fn, stDict, *args)
        Adds a reaction containing a rate function and stoichiometry, that depends on additional information (time).
    rates(y, **kwargs)
        Returns the model rates at time point zero if not specified else in kwargs
    ratesArray(y, **kwargs)
        Returns the model rates at time point zero if not specified else in kwargs
    add_algebraicModule(convert_func, module_name, cpds, derived_cpds)
        Adds an algebraic module
    allCpdNames()
        Returns a list of all compound names, including the compounds derived from algebraic modules
    model(y, t, **kwargs)
        Returns right-hand side of the ODE system
    fullConcVec(z)
        Returns the concentration vector, including all derived variables from algebraic modules
    stoichiometryMatrix()
        Returns the stoichiometric matrix
    print_stoichiometryMatrix()
        Returns a pandas DataFrame of the stoichiometric matrix
    print_stoichiometries()
        Prints stoichiometries
    print_stoichiometries_by_compounds()
        Prints stoichiometries ordered by compounds
    """

    @staticmethod
    def idx(list):
        """ Enumerates a list

        Parameters
        ----------
        list : list
            The list to be enumerated

        Returns
        -------
        enumerated list : dict
            A dictionary containing the list items as keys and the index as values.
        """
        return {it: id for id, it in enumerate(list)}

    def __init__(self, pars={}):
        self.par = ParameterSet(pars)
        self.cpdNames = []
        self.rateFn = {}
        self.stoichiometries = {}
        self.cpdIdDict = {}
        self.algebraicModules = {}

    def store(self, filename):
        """Stores the ParameterSet object into a pickle file

        Parameters
        ----------
        filename : str
            The name of the pickle file

        Returns
        -------
        None
        """
        f = open(filename,'wb')
        pickle.dump(self.par, f)
        f.close()

    @classmethod
    def load(cls, filename):
        """ Loads parameters from a pickle file and stores them in a new model instance

        Parameters
        ----------
        filename : str
            The name of the pickle file

        Returns
        -------
        model : modelbase.model
            A model instance containing the loaded parameters
        """
        par = pickle.load(open(filename,'rb'))
        m = cls(pars=par.__dict__)
        return m


    def rateNames(self):
        """Returns a list of the rate names

        Returns
        -------
        rate_names : list
            A list of the rate names
        """
        return list(self.stoichiometries.keys())


    def updateCpdIds(self):
        """Updates the compound ID dictionary.
        Only needed after modification of the model structure.

        Returns
        -------
        None
        """
        cpdIdDict = self.idx(self.cpdNames)
        cnt = len(self.cpdNames)

        for ammod in self.algebraicModules.values():
            cpdIdDict.update({it: id for id, it in enumerate(ammod['derived_cpds'], cnt)})
            cnt += len(ammod['derived_cpds'])

        self.cpdIdDict = cpdIdDict

    def set_cpds(self, cpdList):
        """Sets the compounds to the given list of compounds.
        This replaces any previously added compounds.

        Parameters
        ----------
        cpdList : list
            A list of compound names (str)

        Returns
        -------
        None
        """
        if isinstance(cpdList, list):
            # Check if all are strings
            if all(isinstance(item, str) for item in cpdList):
                # Check for duplicates
                if len(cpdList) == len(list(set(cpdList))):
                    self.cpdNames = cpdList
                    self.updateCpdIds()
                else:
                    for i in set(cpdList):
                        cpdList.remove(i)
                    raise ValueError("Duplicate entries {}".format(cpdList))
            else:
                raise TypeError("All ist entries have to be strings")
        else:
            raise TypeError("Function expects list input")

    def add_cpd(self, cpdName):
        """Adds a compound to the model.

        Parameters
        ----------
        cpdName : str
            Name of the compound

        Returns
        -------
        None
        """
        if isinstance(cpdName, str):
            if cpdName not in self.cpdIds():
                self.cpdNames.append(cpdName)
                self.updateCpdIds()
            else:
                raise ValueError("Duplicate entry {}".format(cpdName))
        else:
            raise TypeError("Function expects string input")

    def add_cpds(self, cpdList):
        """Adds a list of compounds to the model

        Parameters
        ----------
        cpdList : list
            List of compounds (str) to be added

        Returns
        -------
        None
        """
        if isinstance(cpdList, list):
            for k in cpdList:
                self.add_cpd(k)
        else:
            raise TypeError("Function expects list input")

    def stoichiometryMatrix(self):
        """Returns the stoichiometric matrix.

        Returns
        -------
        stoichiometric matrix : numpy.matrix
        """
        cid = self.idx(self.cpdNames)
        rn = self.rateNames()
        N = np.zeros([len(self.cpdNames),len(rn)])
        for i in range(len(rn)):
            for (c, n) in self.stoichiometries[rn[i]].items():
                N[cid[c],i] = n
        return np.matrix(N)

    def print_stoichiometryMatrix(self):
        """Returns a pandas DataFrame of the stoichiometric matrix.

        Returns
        -------
        stoichiometric matrix : pandas.DataFrame
        """
        M = self.stoichiometryMatrix()
        return pd.DataFrame(M, self.cpdNames, self.rateNames())

    def cpdIds(self):
        """Returns a dictionary with the compound names as keys and their index as values.

        Returns
        -------
        cpdIdDict : dict
        """
        return self.cpdIdDict


    def get_argids(self, *args):
        """Returns the compound IDs for the given args.

        Parameters
        ----------
        *args : list
            List containing the compound names (str)

        Returns
        -------
        cpd_ids : numpy.array
            Numpy array of compound IDs.
        """
        cids = self.cpdIds()
        return np.array([cids[x] for x in args])

    def find_re_argids(self, regexp):
        """Returns indices for compound names matching the regular expression.

        Useful especially in conjunction with labelModel:
        e.g. find all FBPs labelled at pos 3: find_re_argids("\AFBP...1..\Z")

        Parameters
        ----------
        regexp : str
            Regular expression for compound names

        Returns
        -------
        cpd_ids : numpy.array
            Numpy array of compound IDs.
        """
        cids = self.cpdIds()
        reids = []
        for cpdName in cids.keys():
            if re.match(regexp,cpdName):
                reids.append(cids[cpdName])
        return np.array(reids)


    def set_rate(self, rateName, fn, *args):
        """Sets a rate.

        Parameters
        ----------
        rateName : str
            Name of the rate
        fn : method
            Rate function
        *args : str
            Compounds of the function. Must be present in model.

        Returns
        -------
        None
        """
        if not isinstance(rateName, str):
            raise TypeError("RateName must be string")

        sids = self.get_argids(*args)
        if len(sids) == 0:
            def v(y,**kwargs):
                return fn(self.par)
        else:
            def v(y,**kwargs):
                cpdarg = y[sids]
                return fn(self.par,*cpdarg)
        self.rateFn[rateName] = v

    def set_ratev(self, rateName, fn, *args):
        """Sets a rate that depends on additional information (time).

        Difference to set_rate: the rate is called with an additional variable **kwargs,
        that contains time as the key 't'.

        Parameters
        ----------
        rateName : str
            Name of the rate
        fn : method
            Rate function
        *args : str
            Variables

        Returns
        -------
        None
        """
        if not isinstance(rateName, str):
            raise TypeError("RateName must be string")

        sids = self.get_argids(*args)
        if len(sids) == 0:
            def v(y,**kwargs):
                return fn(self.par,**kwargs)
        else:
            def v(y,**kwargs):
                cpdarg = y[sids]
                return fn(self.par,*cpdarg,**kwargs)
        self.rateFn[rateName] = v


    def set_stoichiometry(self, rateName, stDict):
        """Sets the reaction stoichiometry for a rate.

        Parameters
        ----------
        rateName : str
            Name of the rate
        stDict : dict
            Dictionary containing the reaction stoichiometry

        Returns
        -------
        None
        """
        if not isinstance(rateName, str):
            raise TypeError("RateName must be string")

        if not isinstance(stDict, dict):
            raise TypeError("stDict must be dictionary")
        self.stoichiometries[rateName] = stDict

    def set_stoichiometry_byCpd(self, cpdName, stDict):
        """Sets the reaction stoichiometry of a rate by compound name.

        Parameters
        ----------
        cpdName : str
            Name of the compound
        stDict : dict
            Dictionary containing the reaction stoichiometry

        Returns
        -------
        None
        """
        if not isinstance(cpdName, str):
            raise TypeError("cpdName must be string")

        if not isinstance(stDict, dict):
            raise TypeError("stDict must be dictionary")

        for k,v in stDict.items():
            if k not in self.stoichiometries:
                self.stoichiometries[k] = {}
            self.stoichiometries[k][cpdName] = v

    def add_reaction(self,rateName, fn, stDict,*args):
        """Adds a reaction containing a rate function and stoichiometry.

        Parameters
        ----------
        rateName : str
            Name of the rate
        fn : method
            Rate function
        stDict : dict
            Dictionary containing the reaction stoichiometry

        Returns
        -------
        None
        """
        self.set_rate(rateName, fn, *args)
        self.set_stoichiometry(rateName, stDict)

    def add_reaction_v(self,rateName, fn, stDict,*args):
        """Adds a reaction containing a rate function and stoichiometry,
        that depends on additional information (time).

        Difference to add_reaction: the rate is called with an additional variable **kwargs,
        that contains time as the key 't'.

        Parameters
        ----------
        rateName : str
            Name of the rate
        fn : method
            Rate function
        stDict : dict
            Dictionary containing the reaction stoichiometry

        Returns
        -------
        None
        """
        self.set_ratev(rateName, fn, *args)
        self.set_stoichiometry(rateName, stDict)

    def rates(self, y, **kwargs):
        """Returns the model rates at time point zero if not specified else in kwargs.

        Parameters
        ----------
        y : numpy.array
            Array containing all compound concentrations
        t : int or float, optional
            kwargs of the simulation time

        Returns
        -------
        rates : dict
            Dictionary of rate names as keys and rates as values
        """
        z = self.fullConcVec(y)
        return {r:self.rateFn[r](z, **kwargs) for r in self.stoichiometries.keys()}

    def ratesArray(self, y, **kwargs):
        """Returns the model rates at time point zero if not specified else in kwargs.

        Parameters
        ----------
        y : numpy.array
            Array containing all compound concentrations
        t : int or float, optional
            kwargs of the simulation time

        Returns
        -------
        rates : numpy.array
            Numpy array of rates. Ordered like stoichiometries dictionary
        """
        v = self.rates(y, **kwargs)
        return np.array([v[k] for k in self.stoichiometries.keys()])

    def model(self, y, t, **kwargs):
        """Returns right-hand side of the ODE system

        Parameters
        ----------
        y : numpy.array
            Compound concentrations
        t : int or float
            Time point
        **kwargs : dict
            User defined functions for dynamic rate function dependencies

        Returns
        -------
        dydt : numpy.array
            Array of all temporal changes required for ODE integration
        """
        dydt = np.zeros(len(y))
        kwargs.update({'t':t})
        v = self.rates(y, **kwargs)
        idx = self.cpdIds()
        for rate,st in self.stoichiometries.items():
            for cpd,n in st.items():
                dydt[idx[cpd]] += n * v[rate]
        return dydt



    def fullConcVec(self, z):
        """Returns the concentration vector, including all derived variables from algebraic modules.

        Parameters
        ----------
        z : numpy.array or list
            Array of compound concentrations

        Returns
        -------
        z : numpy.array
            Array of compound concentrations and derived variables
        """
        for ammod in self.algebraicModules.values():
            zam = ammod['convert_func'](z)
            z = np.hstack([z,zam])
        return z

    def print_stoichiometries(self):
        """Prints stoichiometries

        Returns
        -------
        None
        """
        pprint.pprint(self.stoichiometries)

    def print_stoichiometries_by_compounds(self):
        """Prints stoichiometries ordered by compounds

        Returns
        -------
        None
        """
        flipped = defaultdict(dict)
        for key, val in self.stoichiometries.items():
            for subkey, subval in val.items():
                flipped[subkey][key] = subval
        pprint.pprint(dict(flipped))

    def add_algebraicModule(self, convert_func, module_name, cpds, derived_cpds):
        """Adds an algebraic module.

        Parameters
        ----------
        module_name : str
            Name of the algebraic module
        convert_func : method
            Method calculating derived variables
        cpds : list
            List of substrates
        derived_cpds : list
            List of compounds derived from the module

        Returns
        -------
        None
        """
        sids = self.get_argids(*cpds)
        def _amwrapper(y):
            if len(np.array(y).shape) == 1:
                cpdarg = y[sids]
                return convert_func(self.par,cpdarg)
            else:
                cpdarg = y[:,sids]
                return np.array([convert_func(self.par,cpdarg[i,:]) for i in range(cpdarg.shape[0])])
        self.algebraicModules[module_name] = {
                'convert_func': _amwrapper,
                'cpds': cpds,
                'derived_cpds': derived_cpds
                }
        self.updateCpdIds()

    def allCpdNames(self):
        """Returns a list of all compound names, including the compounds
        derived from algebraic modules

        Returns
        -------
        names : list
            List of compound names
        """
        names = []
        names.extend(self.cpdNames)
        for ammod in self.algebraicModules.values():
            names.extend(ammod['derived_cpds'])
        return names


class LabelModel(Model):
    """Extension of model class able to create models with
    carbon labeling pattern information.

    Attributes
    ----------
    cpdBaseNames : dict
        Dictionary containing the compounds as keys and number of carbons as values
    Methods
    -------
    add_base_cpd(cpdName, c)
        Adds a carbon containing compound to the model
    add_carbonmap_reaction(rateBaseName, fn, carbonmap, sublist, prodList, *args, **kwargs)
        Sets rates for all reactions of all isotope labelling patterns of the substrates
    set_initconc_cpd_labelpos(y0dict, labelpos)
        Generates a vector of inital concentrations for all label patterns
    """

    @staticmethod
    def generateLabelCpds(cpdName, c):
        """Generates label versions of a compound.

        Adds a string of ones (unlabeled) and zeros (labeled)
        to the compound name. So a compound ABC with two carbons
        can have 4 (in general 2**n) possible patterns:
        ABC00
        ABC01
        ABC10
        ABC11

        Parameters
        ----------
        cpdName : str
            Name of the compound
        c : int
            Number of carbon atoms of the compound

        Returns
        -------
        cpdList : list
            List of compound names with all labeling patterns
        """
        cpdList = [cpdName+''.join(i) for i in itertools.product(('0','1'), repeat = c)]
        return cpdList

    @staticmethod
    def mapCarbons(sublabels, carbonmap):
        """Generates a redistributed string for the substrates, according to the carbonmap

        Parameters
        ----------
        sublabels : list
            List of substrate label positions
        carbonmap : list
            Carbon transition map

        Returns
        -------
        prodLabels : str
            Redistributed string for the substrates (sublabels) according to carbonmap
        """
        prodlabels = ''.join([sublabels[carbonmap[i]] for i in range(len(carbonmap))])
        return prodlabels

    @staticmethod
    def splitLabel(label, numc):
        """Splits the label string according to the lengths given in numc

        Parameters
        ----------
        numc : int
            Number of carbons

        Returns
        -------
        splitLabels : list
            List of label positions
        """
        splitlabels = []
        cnt = 0
        for i in range(len(numc)):
            splitlabels.append(label[cnt:cnt+numc[i]])
            cnt += numc[i]
        return splitlabels

    def __init__(self, pars={}):
        super(LabelModel,self).__init__(pars)
        self.cpdBaseNames = {}

    def add_base_cpd(self, cpdName, c):
        """Adds a carbon containing compound to the model

        Parameters
        ----------
        cpdName : str
            Name of the compound
        c : int
            Number of carbon atoms of the compound

        Returns
        -------
        None
        """
        self.cpdBaseNames[cpdName] = c
        labelNames = self.generateLabelCpds(cpdName,c)
        super(LabelModel,self).add_cpds(labelNames) # add all labelled names

        # now define an algebraic module for the sum of all labels
        # e.g. if CO20, CO21 are the unlabelled and labelled CO2's,
        # the total can be accessed by 'CO2' (likewise for any other more complicated compound)
        if c > 0:
            cpdTotalName = cpdName + '_total'
            def totalconc(par, y):
                return np.array([y.sum()])
            self.add_algebraicModule(totalconc,cpdTotalName, labelNames,[cpdTotalName])


    def add_carbonmap_reaction(self, rateBaseName, fn, carbonmap, subList, prodList, *args, **kwargs):
        """Sets rates for all reactions of all isotope labelling patterns of the substrates.

        Parameters
        ----------
        rateBaseName : str
            Name of the rate
        fn : method
            Rate function
        carbonmap : list
            Carbon transition map
        subList : list
            List of the substrates
        prodList : list
            List of the products
        args : str
            Function arguments (e.g. substrates)
        extLabels : str, optional
            kwargs to define which labels should be added if the number of
            carbon atoms of the products exceeds the number of carbon atoms
            of the substrates. Defaults to all positions labeled.

        Returns
        -------
        None
        """
        # first collect the lengths (num of C) of the substrates and products
        cs = np.array([self.cpdBaseNames[s] for s in subList])
        cp = np.array([self.cpdBaseNames[p] for p in prodList])
        # get all args from *args that are not substrates (can be passed directly)
        otherargs = list(args[len(cs):len(args)])
        # get all possible combinations of label patterns for substrates
        rateLabels = self.generateLabelCpds('',cs.sum())

        extLabels = ''
        if cp.sum() > cs.sum(): # this means labels are introduced to the system
            nrExtLabels = cp.sum() - cs.sum()
            if 'extLabels' in kwargs:
                extLabelList = ['0'] * nrExtLabels
                for extL in kwargs['extLabels']:
                    extLabelList[extL] = '1'
                extLabels = ''.join(extLabelList)
            else:
                extLabels = '1' * (cp.sum() - cs.sum())

        for l in rateLabels: # loop through all patterns
            pl = self.mapCarbons(l+extLabels, carbonmap) # get product labels
            sublabels = self.splitLabel(l, cs)
            prodlabels = self.splitLabel(pl, cp)

            subargs = [args[i]+sublabels[i] for i in range(len(cs))]
            prodargs = [prodList[i]+prodlabels[i] for i in range(len(cp))]
            rateName = rateBaseName+l
            # set rate
            rateargs = subargs+otherargs
            self.set_rate(rateName, fn, *rateargs)
            # set stoichiometry dictionary
            # FIXME think about the possibility that a stoichiometry is not +/-1...
            stDict = {k:-1 for k in subargs}
            for k in prodargs:
                if k in stDict:
                    stDict[k] += 1
                else:
                    stDict[k] = 1
            self.set_stoichiometry(rateName, stDict)

    def set_initconc_cpd_labelpos(self, y0dict, labelpos={}):
        """Generates a vector of inital concentrations for all label patterns.
        Defaults to unlabeled compounds, except those specified in labelpos.

        Parameters
        ----------
        y0dict : dict
            Dictionary with compound names as keys and concentrations as values
        labelpos : dict
            Dictionary with compound names as keys and position of label as value

        Returns
        -------
        y0 : numpy.array
            Initial concentration array for all label patterns
        """
        y0 = np.zeros(len(self.cpdNames))
        for cpd, c in self.cpdBaseNames.items():
            labels = ['0'] * c
            if cpd in labelpos:
                labels[labelpos[cpd]] = '1'
            cpdName = cpd+''.join(labels)
            y0[self.get_argids(cpdName)] = y0dict[cpd]

        return y0
