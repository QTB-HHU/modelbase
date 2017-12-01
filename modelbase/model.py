__author__ = 'oliver'



import modelbase.parameters
from modelbase.algebraicModule import AlgebraicModule

import numpy as np

import scipy.optimize as opt

import numdifftools as nd

import itertools

import re

import pickle



class Model(object):
    '''The base class for modelling. Provides basic functionality.

    This class defines an object with which model construction and
    numeric simulations are made easy.

    An instance of class Model is used to define the model, i.e. the
    dynamic variables and the dynamic equations defining their temporal
    derivatives.

    The numeric simulation is performed with an instance of class
    Simulate.

    Useful analysis methods are provided by class Results

    Mini tutorial
    =============

    Every model is defined by 
    - model parameters
    - model variables
    - rate equations
    - stoichiometries

    Example: A chemical reaction chain 

    -> X -> Y ->
    
    Two variables "X", "Y"

    Three parameters: influx (v0), rate constant conversion X->Y (k1),
    rate constant for outflux (k2)

    Three rate equations:
    - v0 (constant)
    - v1 = k1*X
    - v2 = k2*Y

    with the stoichiometries
    - v0 adds one X
    - v1 removes one X, adds one Y
    - v2 removes one Y

    Mathematically, this results in the two model equations:
    - dX/dt = v0 - k1*X
    - dY/dt = k1*X - k2*Y

    When instanciating a model, the model parameters are provided as
    a dictionary:

    m = Model({'v0':1, 'k1': 0.5, 'k2': 0.1})

    The variables can now be accessed by m.par.v0, m.par.k1 and m.par.k2

    Now, the variables need to be added. Variables are ALWAYS defined by
    names (i.e. strings). These are later used to access and identify
    the variables and their values. Here:

    m.set_cpds(['X','Y'])

    The last thing to do is to set the rates. This is done using
    set_rate.  Here, the first argument is always a name that the rate
    is associated with (to access it later) and the second is a
    function that calculates the rate. The remaining parameters are
    the names of the variables whose values are passed to the
    function. The function must always accept as first argument a
    parameter object (actually, m.par), and the remaining arguments
    are the values of the variables used to calculate the rate.

    Rate v0: this is particularly simple, because it is constant:

    m.set_rate('v0', lambda p: p.v0)

    Rate v1: this depends also on the value of variable 'X'. So we define
    a function first.

    def v1(p,x):
        return p.k1*x

    m.set_rate('v1', v1, 'X')

    Likewise v2:

    def v2(p,x):
        return p.k2*x

    m.set_rate('v2', v2, 'Y')

    Last thing is to set the stoichiometries:

    m.set_stoichiometry('v0',{'X':1})
    
    m.set_stoichiometry('v1',{'X':-1,'Y':1})
    
    m.set_stoichiometry('v2',{'Y':-1})

    Simulation and Plot:

    s = Simulate(m)
    
    T = np.linspace(0,100,1000)
    Y = s.timeCourse(T,np.zeros(3))

    plt.plot(T,Y)

    This example is found in example.py. Other examples using additional
    functionalities are provided in the other example{i}.py files
    '''

    @staticmethod
    def idx(list):
        return {it: id for id, it in enumerate(list)}


    def __init__(self, pars={}, defaultpars={}):
        self.par = modelbase.parameters.ParameterSet(pars,defaultpars)
        self.cpdNames = []
        self.rateFn = {}
        self.stoichiometries = {}
        self.cpdIdDict = {}


    def store(self, filename):
        '''
        stores the parameters to file FILENAME
        :input filename: FILENAME
        '''
        f = open(filename,'wb')
        pickle.dump(self.par, f)
        f.close()

    @classmethod
    def load(cls, filename):
        '''
        loads parameters from file and invokes constructor of corresponding class
        :input filename: where to load parameters from
        '''
        par = pickle.load(open(filename,'rb'))
        m = cls(pars=par.__dict__)
        return m


    def rateNames(self):
        return list(self.stoichiometries.keys())


    def updateCpdIds(self):
        '''
        updates self.cpdIdDict. Only needed after modification of model
        structure, e.g. by set_cpds, add_cpd and add_cpds
        '''
        self.cpdIdDict = self.idx(self.cpdNames)

    def set_cpds(self,cpdList):
        '''
        sets the names of the compounds in model to cpdList
        TO DECIDE: do we want add_cpds and even rm_cpds?
        '''
        self.cpdNames = cpdList
        self.updateCpdIds()

    def add_cpd(self, cpdName):
        '''
        adds a single compound with name cpdName (string) to cpdNames
        if it does not yet exist
        '''
        if cpdName not in self.cpdIds():
            self.cpdNames.append(cpdName)
            self.updateCpdIds()

    def add_cpds(self, cpdList):
        '''
        adds a list of compounds (list of strings with names) to cpdNames
        '''
        for k in cpdList:
            self.add_cpd(k)
        #self.cpdNames = self.cpdNames + cpdList
        #self.updateCpdIds()

    def stoichiometryMatrix(self):
        '''
        returns the stoichiometry matrix
        '''

        cid = self.idx(self.cpdNames)
        #print cid
        rn = self.rateNames()
    
        N = np.zeros([len(self.cpdNames),len(rn)])

        for i in range(len(rn)):
            for (c, n) in self.stoichiometries[rn[i]].items():
                #print "c=%s, cid=%d, r=%s, n=%d" % (c, cid[c], rn[i], n)
                N[cid[c],i] = n

        return np.matrix(N)



    def cpdIds(self):
        '''
        returns a dict with keys:cpdNames, values:idx
        This is now cached in self.cpdIdDict. This is updated whenever
        a compound is added.
        '''
        return self.cpdIdDict
        #return {self.cpdNames[i]:i for i in range(len(self.cpdNames))}


    def get_argids(self, *args):

        cids = self.cpdIds()
        return np.array([cids[x] for x in args])

    def find_re_argids(self, regexp):
        '''
        Returns list of indices for which the compound name matches the
        regular expression
        Useful especially in conjunction with labelModel:
        e.g. find all FBPs labelled at pos 3: find_re_argids("\AFBP...1..\Z")
        '''
        cids = self.cpdIds()
        reids = []
        for cpdName in cids.keys():
            if re.match(regexp,cpdName):
                reids.append(cids[cpdName])
        return np.array(reids)


    def set_rate(self, rateName, fn, *args):
        '''
        sets a rate. Arguments:
        Input: rateName (string), fn (the function) and _names_ of compounds which are passed to the function.
        The function fn is called with the parameters self.par as first argument and the dynamic variables corresponding to the compounds as variable argument list.

        Example
        -------
        m = modelbase.model.Model({'k1':0.5})
        m.set_cpds(['X','Y','Z'])
        def v1(par,x):
            return par.k1*x
        m.set_rate('v1',v1,'X')

        m.rateFn['v1'](np.array([3,2,1]))
        # 1.5
        '''
        
        sids = self.get_argids(*args)


        if len(sids) == 0:
            # note: the **kwargs is necessary to allow all rates to be called in the same way. It can be empty.
            def v(y,**kwargs): 
                return fn(self.par)
        else:
            def v(y,**kwargs):
                cpdarg = y[sids]
                return fn(self.par,*cpdarg)

        self.rateFn[rateName] = v



    def set_ratev(self, rateName, fn, *args):
        '''
        sets a rate, which depends on additional information. 
        Difference to set_rate: the rate is called with an additional variable **kwargs. 
        This always contains time as key 't', and other user-defined stuff that is passed to methods 'model', 'rates'
        Arguments:
        Input: rateName (string), fn (the function) and _names_ of compounds which are passed to the function.
        The function fn is called with the parameters self.par as first argument and the dynamic variables corresponding to the compounds as variable argument list.

        Example
        -------
        m = modelbase.model.Model({'l':1,'k1':0.5})
        m.set_cpds(['X'])
        def v1(par,**kwargs):
            return np.exp(-par.l*kwargs['t'])
        m.set_ratev('v1',v1)

        m.rateFn['v1'](np.array([0]),t=0)
        # 1
        m.rateFn['v1'](np.array([0]),t=1)
        # 0.36787944117144233
        '''
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
        '''
        sets stoichiometry for rate rateName to values contained in stDict
        
        Example
        -------
        m.set_stoichiometry('v1',{'X':-1,'Y',1})

        '''

        self.stoichiometries[rateName] = stDict


    def set_stoichiometry_byCpd(self, cpdName, stDict):
        '''
        same as set_stoichiometry, but by compound name
        '''
        for k,v in stDict.items():
            self.stoichiometries[k][cpdName] = v



    def rates(self, y, **kwargs):
        '''
        argument: np.array y - values of all compounds
        output: dict with rateNames as keys and corresponding values
        '''

        return {r:self.rateFn[r](y, **kwargs) for r in self.stoichiometries.keys()}


    def ratesArray(self, y, **kwargs):
        '''
        argument: np.array y - values of all compounds
        output: array with rates, order as self.stoichiometry.keys()
        '''

        v = self.rates(y, **kwargs)
        return np.array([v[k] for k in self.stoichiometries.keys()])
        
        

    def model(self, y, t, **kwargs):
        '''
        argument: np.array y - including values of all compounds
        output: np.array dydt - including all corresponding temporal changes required for dynamic simulation / ODE integration
        '''

        dydt = np.zeros(len(y))

        kwargs.update({'t':t})

        v = self.rates(y, **kwargs)
        idx = self.cpdIds()

        for rate,st in self.stoichiometries.items():
            for cpd,n in st.items():
                dydt[idx[cpd]] += n * v[rate]

        return dydt



    def fullConcVec(self, y):
        '''
        included only for compatibility reasons. Now an AlgmSimulate can be used also if no algebraic module is present
        '''
        return y


    def numericElasticities(self, y0, rate):
        '''
        y0: state vector
        rate: name of rate for which elasticities shall be determined
        '''
 
        def vi(y):
            v = self.rates(y)
            return v[rate]

        jac = nd.Jacobian(vi,step=y0.min()/100)

        epsilon = jac(y0)

        return epsilon

    def allElasticities(self, y0, norm=False):
        ''' 
        calculates all elasticities:
        :param y0: state vector
        :return: all elasticities as np.matrix
        '''

        rateIds = self.rateNames()

        epsilon = np.zeros([len(rateIds), len(self.cpdNames)])

        for i in range(len(rateIds)):

            def vi(y):
                return self.rateFn[rateIds[i]](y)
                
            jac = nd.Jacobian(vi, step=y0.min()/100)

            epsilon[i,:] = jac(y0)

        if norm:
            v = np.array(self.rates(y0).values())
            epsilon = (1/v).reshape(len(v),1)*epsilon*y0

        return np.matrix(epsilon)


    def numericJacobian(self, y0, **kwargs):
        '''
        y0: state vector at which Jacobian is calculated
        '''
        J = np.zeros([len(y0),len(y0)])

        if np.isclose(y0.min(),0):
            jstep = None
        else:
            jstep = y0.min()/100

        for i in range(len(y0)):

            def fi(y):
                dydt = self.model(y, 0, **kwargs)
                return dydt[i]

            jac = nd.Jacobian(fi,step=jstep)

            J[i,:] = jac(y0)

        return np.matrix(J)


    def findSteadyState(self, y0, **kwargs):
        '''
        tries to find the steady-state by numerically solving the algebraic system dy/dt = 0.
        input: y0: initial guess
        TODO: this method can be improved. So far, it simply tries the standard solving method hybr
        '''
        
        def fn(x):
            return self.model(x, 0, **kwargs)
        sol = opt.root(fn, y0)

        if sol.success == True:
            return sol.x
        else:
            return False



    def concentrationControlCoefficients(self, y0, pname, norm=True, **kwargs):
        '''
        invokes findSteadyState to calculate the concentration control coefficients
        for parameter pname
        :input y0: initial guess for steady-state
        :input pname: parameter name to vary
        :input norm: if True (default), normalize coefficients
        :returns: response coefficients
        '''

        origValue = getattr(self.par, pname)

        def fn(x):
            self.par.update({pname: x})
            return self.findSteadyState(y0, **kwargs)

        jac = nd.Jacobian(fn, step=origValue/100.)

        cc = np.array(jac(origValue))

        self.par.update({pname: origValue})

        if norm:
            ss = self.findSteadyState(y0, **kwargs)
            cc = origValue * cc / ss.reshape(ss.shape[0],1)

        return cc




###### class AlgmModel #################################################################

class AlgmModel(Model):
    '''
    Subclass of Model, which incorporates algebraic modules.
    An algebraic module is basically a function that allows calculation of concentrations
    of some variables by other variables. 
    The simplest example is a conserved quantity, e.g. ATP+ADP=Atotal, then ADP=Atotal-ATP
    can be determined from ATP.
    Rapid Equilibrium modules are a typical application.
    '''
    def __init__(self, pars={}, defaultpars={}):
        super(AlgmModel,self).__init__(pars,defaultpars)
        self.algebraicModules = []


    def updateCpdIds(self):
        '''
        updates self.cpdIdDict. Only needed after modification of model
        structure, e.g. by set_cpds, add_cpd and add_cpds
        '''
        cpdIdDict = self.idx(self.cpdNames)
        cnt = len(self.cpdNames)

        for ammod in self.algebraicModules:
            cpdIdDict.update({it: id for id, it in enumerate(ammod['amCpds'], cnt)})
            cnt += len(ammod['amCpds'])
        
        self.cpdIdDict = cpdIdDict
        

    def add_algebraicModule(self, am, amVars, amCpds):
        '''
        this adds a module in which several compound concentrations can be calculated algebraicly.
        am: modelbase.algebraicModule.AlgebraicModule
        amVars: list of names of variables used for module in embedding model
        amCpds: list of names of compounds which are calculated by the module from amVars
        '''
        
        self.algebraicModules.append({'am': am, 'amVars': amVars, 'amCpds': amCpds})
        self.updateCpdIds()

    
    #def get_argids(self, *args):
    #    # FIXME: this should also be cached
    #    cpdids = {it: id for id, it in enumerate(self.cpdNames)}
    #    cnt = len(self.cpdNames)
    #
    #    for ammod in self.algebraicModules:
    #        cpdids.update({it: id for id, it in enumerate(ammod['amCpds'], cnt)})
    #        cnt += len(ammod['amCpds'])
    #
    #    return np.array([cpdids[x] for x in args])

    """
    def set_rate(self, rateName, fn, *args):
        '''
        sets a rate. Arguments:
        Input: rateName (string), fn (the function) and _names_ of compounds which are passed to the function.
        The function fn is called with the parameters self.par as first argument and the dynamic variables corresponding to the compounds as variable argument list.

        In contrast to class Model, args can contain a name from an algebraic module

        '''

        cpdids = {it: id for id, it in enumerate(self.cpdNames)}
        cnt = len(self.cpdNames)

        for ammod in self.algebraicModules:
            cpdids.update({it: id for id, it in enumerate(ammod['amCpds'], cnt)})
            cnt += len(ammod['amCpds'])

        argids = np.array([cpdids[x] for x in args])

        if len(argids) == 0:
            def v(y):
                return fn(self.par)
        else:
            def v(y):
                cpdarg = y[argids]
                return fn(self.par,*cpdarg)

        self.rateFn[rateName] = v
    """    

    def fullConcVec(self, y):
        '''
        returns the full concentration vector, including all concentrations from algebraic modules
        input: y - state vector of all dynamic variables
        output: z - state vector extended by all derived concentrations
        '''
        z = y.copy()
        #vlist = [y]
        #cpdids = {it: id for id, it in enumerate(self.cpdNames)}
        cpdids = self.cpdIds()

        for ammod in self.algebraicModules:
            varids = np.array([cpdids[x] for x in ammod['amVars']])
            if len(z.shape) == 1:
                zin = z[varids]
            else:
                zin = z[:,varids]
            zam = ammod['am'].getConcentrations(zin)

            #cpdidsam = {it:id for id,it in enumerate(ammod['amCpds'], z.size)}

            #if len(zam.shape) == 1:
            #    zam = zam[:,np.newaxis]

            z = np.hstack([z,zam])
            #cpdids = dict(cpdids, **cpdidsam)
            #vlist.append(ammod['am'].getConcentrations(y[varids]))

        #z = np.hstack(vlist)

        return z


    def rates(self, y, **kwargs):

        z = self.fullConcVec(y)

        return {r:self.rateFn[r](z, **kwargs) for r in self.stoichiometries.keys()}



    def allCpdNames(self):
        ''' returns list of all compounds, including from algebraic modules '''
        names = []
        names.extend(self.cpdNames)
        for ammod in self.algebraicModules:
            names.extend(ammod['amCpds'])

        return names

 
    def allElasticities(self, y0, norm=False):
        ''' 
        calculates all _direct_ elasticities:
        Rates usually depend on a concentration and not directly on a conserved equilbrium module variable.
        Therefore, the partial derivatives of the rate expression itself is zero wrt the equilibrium variable, but non-zero wrt to the concentration.
        :param y0: state vector
        :return: all elasticities as np.matrix
        # FIXME: more elegant merge with superclass method
        '''
    
        rateIds = self.rateNames()
    
        epsilon = np.zeros([len(rateIds), len(self.allCpdNames())])
    
        z0 = self.fullConcVec(y0)
    
        for i in range(len(rateIds)):
    
            def vi(y):
                return self.rateFn[rateIds[i]](y)
                
            jac = nd.Jacobian(vi, step=z0.min()/100)
    
            epsilon[i,:] = jac(z0)
    
        if norm:
            v = np.array(self.rates(z0).values())
            epsilon = (1/v).reshape(len(v),1)*epsilon*z0

        return np.matrix(epsilon)







###### class LabelModel #################################################################


#def generateLabelCpds(cpdName, c):
#    '''
#    generates label versions of a compound.
#    input: string cpdName, int c (number of carbon atoms)
#    output: list of compounds names with all labeling patterns accroding to Name000, Name001 etc
#    '''
#
#    cpdList = [cpdName+''.join(i) for i in itertools.product(('0','1'), repeat = c)]
#
#    return cpdList
#
#def mapCarbons(sublabels, carbonmap):
#    '''
#    generates a redistributed string for the substrates (sublabels) according to carbonmap
#    '''
#    prodlabels = ''.join([sublabels[carbonmap[i]] for i in range(len(carbonmap))])
#    return prodlabels
#
#def splitLabel(label, numc):
#    '''
#    splits the label string according to the lengths given in the list/vector numc
#    '''
#    splitlabels = []
#    cnt = 0
#    for i in range(len(numc)):
#        splitlabels.append(label[cnt:cnt+numc[i]])
#        cnt += numc[i]
#    return splitlabels


class LabelModel(AlgmModel):
    '''
    LabelModel allows to define a model with carbon labelling pattern information
    
    Important information on usage:
    -------------------------------
    Compounds must be added with the add_base_cpd method, which here takes two arguments:
    cpdName (string) and c (int) specifying number of carbon atoms
    '''


    # some important static methods

    @staticmethod
    def generateLabelCpds(cpdName, c):
        '''
        generates label versions of a compound.
        input: string cpdName, int c (number of carbon atoms)
        output: list of compounds names with all labeling patterns accroding to Name000, Name001 etc
        '''
    
        cpdList = [cpdName+''.join(i) for i in itertools.product(('0','1'), repeat = c)]
    
        return cpdList

    @staticmethod    
    def mapCarbons(sublabels, carbonmap):
        '''
        generates a redistributed string for the substrates (sublabels) according to carbonmap
        '''
        prodlabels = ''.join([sublabels[carbonmap[i]] for i in range(len(carbonmap))])
        return prodlabels
    
    @staticmethod
    def splitLabel(label, numc):
        '''
        splits the label string according to the lengths given in the list/vector numc
        '''
        splitlabels = []
        cnt = 0
        for i in range(len(numc)):
            splitlabels.append(label[cnt:cnt+numc[i]])
            cnt += numc[i]
        return splitlabels





    def __init__(self, pars={}, defaultpars={}):
        super(LabelModel,self).__init__(pars,defaultpars)
        self.cpdBaseNames = {}


    def add_base_cpd(self, cpdName, c):
        '''
        adds compound to model, generating all possible labelling patterns
        :param cpdName: compound base name
        :param c: number of C atoms
        '''
        self.cpdBaseNames[cpdName] = c
        labelNames = self.generateLabelCpds(cpdName,c)
        super(LabelModel,self).add_cpds(labelNames) # add all labelled names

        # now define an algebraic module for the sum of all labels
        # e.g. if CO20, CO21 are the unlabelled and labelled CO2's,
        # the total can be accessed by 'CO2' (likewise for any other more complicated compound)
        def totalconc(par, y):
            return np.array([y.sum()])
        tc = AlgebraicModule({}, totalconc)
        self.add_algebraicModule(tc,labelNames,[cpdName])


    def add_carbonmap_reaction(self, rateBaseName, fn, carbonmap, subList, prodList, *args, **kwargs):
        '''
        sets all rates for reactions for all isotope labelling patterns of the substrates.
        Sets all stoichiometries for these reactions.
        requires additionally
        - carbonmap: a list defining how the carbons appear in the products
          (of course, number of Cs must be the same for substrates and products,
           _except_ if 
             1. a uniform outflux is defined. Then simply carbonmap=[])
             2. extra labels enter the system. Then the pattern is defined by kwargs['extLabels']. Defaults to all labelled.
        - subList: list of substrates
        - prodList: list of products
        - *args: list of arguments required to calculate rate using function fn 
          (including substrates and possibly allosteric effectors). 
          In this list, substrate names MUST come first
        - **kwargs: required to
             - define extra lables. Key 'extLabels', Value: list of labels, starting with 0

        examples for carbon maps:
        TPI: GAP [0,1,2] -> DHAP [2,1,0] (order changes here), carbonmap = [2,1,0]
        Ald: DHAP [0,1,2] + GAP [3,4,5] -> FBP, carbonmap = [0,1,2,3,4,5]
        TK: E4P [0,1,2,3] + X5P [4,5,6,7,8] -> GAP [6,7,8] + F6P [4,5,0,1,2,3], carbonmap = [6,7,8,4,5,0,1,2,3]
        '''

        # first collect the lengths (num of C) of the substrates and products
        cs = np.array([self.cpdBaseNames[s] for s in subList])
        cp = np.array([self.cpdBaseNames[p] for p in prodList])

        # get all args from *args that are not substrates (can be passed directly)
        otherargs = list(args[len(cs):len(args)])
        #print "otherargs:", otherargs

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
                extLabels = '1' * (cp.sum() - cs.sum()) # FIXME make more flexible to allow labels and no-labels to be introduced

        for l in rateLabels: # loop through all patterns
            #print l
            pl = self.mapCarbons(l+extLabels, carbonmap) # get product labels
            sublabels = self.splitLabel(l, cs)
            prodlabels = self.splitLabel(pl, cp)

            subargs = [args[i]+sublabels[i] for i in range(len(cs))]
            #print subargs
            prodargs = [prodList[i]+prodlabels[i] for i in range(len(cp))]
            #print prodargs

            rateName = rateBaseName+l

            # set rate
            rateargs = subargs+otherargs
            #print rateargs
            self.set_rate(rateName, fn, *rateargs)

            # set stoichiometry dictionary
            # FIXME think about the possibility that a stoichiometry is not +/-1...
            stDict = {k:-1 for k in subargs}
            for k in prodargs:
                if k in stDict:
                    stDict[k] += 1
                else:
                    stDict[k] = 1
            #stDict.update({k:1 for k in prodargs}) # did not work if substrates = products
            #print stDict
            self.set_stoichiometry(rateName, stDict)
            

    def set_initconc_cpd_labelpos(self, y0dict, labelpos={}):
        '''
        generates a vector of initial concentrations, such that
        everything is unlabelled excpet those specified in dictionary labelpos.
        :param y0dict: dict with compound names as keys and total concentrations as values
        :param labelpos: dict with compound names as keys and the position of the label as value
        :return: A full length vector of concentrations.

        Inputs:

        y0dict: a dictionary with compound names as keys and concentrations as values. These are used to set the total concentrations. By default to the unlabelled compound.

        labelpos: a dictionary with compound names as keys and the position of the label as value. 

        Output:

        A full length vector of concentrations.

        Example: GAP labelled at 1-position, DHAP and FBP unlabelled
        
        y0 = m.set_initconc_cpd_labelpos({'GAP':1,'DHAP':20,'FBP':4},{'GAP':0})
        '''
        y0 = np.zeros(len(self.cpdNames))
        for cpd, c in self.cpdBaseNames.items():
            labels = ['0'] * c
            if cpd in labelpos:
                labels[labelpos[cpd]] = '1'
            cpdName = cpd+''.join(labels)
            y0[self.get_argids(cpdName)] = y0dict[cpd]

        return y0

