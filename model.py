__author__ = 'oliver'


import modelbase.parameters
import numpy as np

import scipy.optimize as opt

import numdifftools as nd


def idx(list):
    return {it: id for id, it in enumerate(list)}


class Model(object):
    '''
    base class for modelling. Provides basic functionality.
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
        '''
        return {self.cpdNames[i]:i for i in range(len(self.cpdNames))}
        

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
        cids = self.cpdIds()
        
        sids = np.array([cids[x] for x in args])

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


    def numericJacobian(self, y0):
        '''
        y0: state vector at which Jacobian is calculated
        '''
        J = np.zeros([len(y0),len(y0)])

        for i in range(len(y0)):

            def fi(y):
                dydt = self.model(y)
                return dydt[i]

            jac = nd.Jacobian(fi,step=y0.min()/100)

            J[i,:] = jac(y0)

        return np.matrix(J)


    def findSteadyState(self, y0):
        '''
        tries to find the steady-state by numerically solving the algebraic system dy/dt = 0.
        input: y0: initial guess
        TODO: this method can be improved. So far, it simply tries the standard solving method hybr
        '''
        
        sol = opt.root(self.model, y0)

        if sol.success == True:
            return sol.x
        else:
            return False
