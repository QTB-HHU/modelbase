__author__ = 'oliver'




import modelbase.parameters
import numpy as np

import scipy.optimize as opt

import numdifftools as nd

import re


def idx(list):
    return {it: id for id, it in enumerate(list)}


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

    def __init__(self, pars={}, defaultpars={}):
        self.par = modelbase.parameters.ParameterSet(pars,defaultpars)
        self.cpdNames = []
        self.rateFn = {}
        self.stoichiometries = {}


    def rateNames(self):
        return self.stoichiometries.keys()


    def set_cpds(self,cpdList):
        '''
        sets the names of the compounds in model to cpdList
        TO DECIDE: do we want add_cpds and even rm_cpds?
        '''
        self.cpdNames = cpdList

    def add_cpd(self, cpdName):
        '''
        adds a single compound with name cpdName (string) to cpdNames
        '''
        self.cpdNames.append(cpdName)

    def add_cpds(self, cpdList):
        '''
        adds a list of compounds (list of strings with names) to cpdNames
        '''
        self.cpdNames = self.cpdNames + cpdList

    def stoichiometryMatrix(self):
        '''
        returns the stoichiometry matrix
        '''

        cid = idx(self.cpdNames)
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
        FIXME: this should be cached to improve efficiency!!!
        '''
        return {self.cpdNames[i]:i for i in range(len(self.cpdNames))}


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
        for cpdName in self.cpdNames:
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
        m.set_rate('v1',v1)

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


    def numericJacobian(self, y0, **kwargs):
        '''
        y0: state vector at which Jacobian is calculated
        '''
        J = np.zeros([len(y0),len(y0)])

        for i in range(len(y0)):

            def fi(y):
                dydt = self.model(y, 0, **kwargs)
                return dydt[i]

            jac = nd.Jacobian(fi,step=y0.min()/100)

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
